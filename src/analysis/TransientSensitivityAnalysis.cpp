/**
 * @author Gabriel Henrique Silva
 * @file TransientSensitivityAnalysis.cpp
 * @brief Transient sensitivity by the direct method, step by step.
 *
 * The analysis hooks into the transient engine through its step observer, so
 * the trajectory is EXACTLY the one a plain `.TRAN` run produces (breakpoints,
 * step rejection and the LTE controller included).
 *
 * Why the DIRECT method here
 * --------------------------
 * The transient sensitivity is the derivative of the whole TRAJECTORY with
 * respect to a parameter, and the trajectory is a recursion:
 *
 *      x_{n+1} = Φ(x_n, p)
 *
 * whose sensitivity obeys
 *
 *      Y_{n+1} = ∂Φ/∂x · Y_n + ∂Φ/∂p
 *
 * The objective at step n is O(x_n), so dO/dp = (∂O/∂x)·Y_n and the
 * sensitivity VECTOR Y_n has to be propagated — one solve per parameter per
 * step. (The adjoint of a recursion is a BACKWARD sweep, which cannot
 * produce the sensitivity of every time point of a forward simulation; the
 * adjoint is the right tool for a single objective at the END of the
 * trajectory, and the direct method is the right tool for a time series.)
 *
 * The recursion is exact for the discrete scheme the simulator integrates:
 *
 *      J_n·Y_n = w_n
 *
 * with J_n = dF/dx + s·dQ/dx (the effective Jacobian of the step) and
 *
 *      w_n = ∂B/∂p − (∂G/∂p)·X_n − s·(∂C/∂p)·X_n     (the device's own
 *            + s·(dQ/dx)·Y_{n-1}                        perturbation)
 *            + s·dQ/dp|_explicit
 *            + e_{n-1}                                 (trapezoidal rule)
 *
 * where the history terms are what make the perturbation PROPAGATE through
 * the trajectory instead of dying at every step:
 *
 *   - s·(dQ/dx)·Y_{n-1} is the TRAJECTORY dependence of the accepted charge
 *     (the same for every device: every device owns disjoint rows of dQ/dx,
 *     so the sum over the devices is one sparse product),
 *   - s·dQ/dp|_explicit is the device's own explicit dependence of its
 *     stored charge on the parameter (reported by the device layer),
 *   - e_{n-1} is the derivative of the trapezoidal rule's secondary history
 *     (the stored dQ/dt), which the analysis carries forward.
 *
 * Cost: one factorization per step (shared by all the parameters) plus one
 * triangular solve per parameter per step.
 */

#include "TransientSensitivityAnalysis.h"

#include "../core/DaeSolution.h"
#include "../core/Device.h"
#include "../core/Node.h"
#include "../devices/Resistor.h"
#include "../solver/ComplexSparseLU.h"
#include "../solver/DaeMatrix.h"
#include "../solver/StateStore.h"
#include "DaeSystem.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace e_XSpice {

namespace {
bool iequals(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i]))) return false;
    }
    return true;
}
} // namespace

TransientSensitivityAnalysis::TransientSensitivityAnalysis(Circuit& circuit)
    : circuit_(circuit) {}

TransientSensitivityAnalysis::~TransientSensitivityAnalysis() = default;

TransientSensitivityResult TransientSensitivityAnalysis::run(bool use_mpi) {
    TransientSensitivityResult result;
    result.output_var = output_var_;

    if (output_var_.empty()) {
        result.message = "no output variable given (use V(node) or I(dev))";
        return result;
    }
    if (tstop_ <= tstart_ || tstep_ <= 0.0) {
        result.message = "invalid time range";
        return result;
    }

    // ---- parameter list ----------------------------------------------------
    std::vector<std::pair<std::string, std::string>> expanded;
    for (const auto& [dev_name, param_name] : params_) {
        Device* dev = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), dev_name)) { dev = dv.get(); break; }
        }
        if (dev && param_name.empty()) {
            for (const auto& pn : dev->sensitivityParameters()) {
                expanded.emplace_back(dev_name, pn);
            }
        } else {
            expanded.emplace_back(dev_name, param_name);
        }
    }
    if (expanded.empty()) {
        for (const auto& d : circuit_.getDevices()) {
            for (const auto& pn : d->sensitivityParameters()) {
                expanded.emplace_back(d->getName(), pn);
            }
        }
    }
    const size_t nparam = expanded.size();
    for (size_t i = 0; i < nparam; ++i) {
        result.param_names.push_back(expanded[i].first + ":" + expanded[i].second);
    }

    // Resolve the devices once (the sweep must not search the device list
    // per parameter per time step).
    std::vector<Device*> param_dev(nparam, nullptr);
    std::vector<bool> param_ok(nparam, false);
    for (size_t i = 0; i < nparam; ++i) {
        Device* dev = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), expanded[i].first)) { dev = dv.get(); break; }
        }
        param_dev[i] = dev;
        if (!dev) continue;
        for (const auto& pn : dev->sensitivityParameters()) {
            if (iequals(pn, expanded[i].second)) { param_ok[i] = true; break; }
        }
    }

    // ---- the DAE layout (node rows first, branch rows after) --------------
    // The DAE numbers the node rows in the circuit's node order (the ground
    // is never materialized) and puts the branch currents after them, so the
    // layout can be reconstructed without asking the engine.
    std::vector<int> node_map;
    int num_node_rows = 0;
    {
        int max_id = 0;
        for (const auto& n : circuit_.getNodes()) {
            max_id = std::max(max_id, n->getId());
        }
        node_map.assign(static_cast<size_t>(max_id) + 1, -1);
        for (const auto& n : circuit_.getNodes()) {
            if (n->isGround()) continue;
            node_map[static_cast<size_t>(n->getId())] = num_node_rows++;
        }
    }
    int num_branch_rows = 0;
    for (const auto& d : circuit_.getDevices()) {
        num_branch_rows += std::max(0, d->numExtraUnknowns());
    }
    const int n_dae = num_node_rows + num_branch_rows;

    // ---- objective gradient (∂O/∂X) in DAE space -------------------------
    std::vector<double> c_dae(static_cast<size_t>(n_dae), 0.0);
    {
        std::string s = output_var_;
        size_t p = s.find_first_not_of(" \t");
        if (p != std::string::npos) s = s.substr(p);
        p = s.find_last_not_of(" \t");
        if (p != std::string::npos) s = s.substr(0, p + 1);
        if (s.size() < 3 || s.back() != ')' || s[1] != '(') {
            result.message = "malformed output variable '" + output_var_ + "'";
            return result;
        }
        const char kind = static_cast<char>(
            std::toupper(static_cast<unsigned char>(s[0])));
        const std::string inner = s.substr(2, s.size() - 3);

        if (kind == 'V') {
            const size_t comma = inner.find(',');
            const std::string na = (comma == std::string::npos) ? inner
                                                              : inner.substr(0, comma);
            const std::string nb = (comma == std::string::npos) ? std::string()
                                                              : inner.substr(comma + 1);
            auto rowOfName = [&](const std::string& name) -> int {
                if (name.empty() || name == "0") return -1;
                for (const auto& nd : circuit_.getNodes()) {
                    if (nd->isGround()) continue;
                    if (!iequals(nd->getName(), name)) continue;
                    const int id = nd->getId();
                    if (id < 0 || id >= static_cast<int>(node_map.size())) return -1;
                    return node_map[static_cast<size_t>(id)];
                }
                return -1;
            };
            const int r1 = rowOfName(na);
            const int r2 = nb.empty() ? -1 : rowOfName(nb);
            if (r1 < 0 && r2 < 0) {
                result.message = "output node of '" + output_var_ + "' not found";
                return result;
            }
            if (r1 >= 0) c_dae[static_cast<size_t>(r1)] += 1.0;
            if (r2 >= 0) c_dae[static_cast<size_t>(r2)] -= 1.0;
        } else {
            Device* d = nullptr;
            for (const auto& dv : circuit_.getDevices()) {
                if (iequals(dv->getName(), inner)) { d = dv.get(); break; }
            }
            if (!d) {
                result.message = "output device '" + inner + "' not found";
                return result;
            }
            if (d->getBranchIndex() >= 0) {
                const int br = d->getBranchIndex();
                if (br >= 0 && br < n_dae) c_dae[static_cast<size_t>(br)] += 1.0;
            } else if (d->getType() == DeviceType::RESISTOR) {
                // I = (v+ - v-)/R  =>  dI/dv+ = 1/R, dI/dv- = -1/R.
                const double gdev =
                    1.0 / static_cast<const Resistor*>(d)->getResistance();
                const Node* np = d->getNodes()[0];
                const Node* nn = d->getNodes()[1];
                if (np && !np->isGround()) {
                    const int id = np->getId();
                    if (id >= 0 && id < static_cast<int>(node_map.size())) {
                        const int r = node_map[static_cast<size_t>(id)];
                        if (r >= 0) c_dae[static_cast<size_t>(r)] += gdev;
                    }
                }
                if (nn && !nn->isGround()) {
                    const int id = nn->getId();
                    if (id >= 0 && id < static_cast<int>(node_map.size())) {
                        const int r = node_map[static_cast<size_t>(id)];
                        if (r >= 0) c_dae[static_cast<size_t>(r)] -= gdev;
                    }
                }
            } else {
                result.message = "I(" + inner + ") is not a supported output";
                return result;
            }
        }
    }

    // ---- CSV header --------------------------------------------------------
    if (csv_) {
        csv_->precision(17);
        if (!csv_header_written_) {
            *csv_ << "time,obj_value";
            for (const auto& name : result.param_names) {
                *csv_ << "," << name;
            }
            *csv_ << "\n";
            csv_header_written_ = true;
        }
    }

    // ---- the adjoint/solve machinery ---------------------------------------
    // The effective matrix changes at every step, but its PATTERN does not.
    // The pivot-safe ordering is decided ONCE (from the first matrix) and
    // frozen; the LU is refilled and refactorized at every step.
    struct Permuted {
        ComplexSparseLU lu;
        std::vector<int> row_of_pos;     // permuted position -> DAE row
        std::vector<int> row_ptr;        // permuted CSR row pointers
        std::vector<int> col_idx;        // permuted CSR column indices
        std::vector<int> src_slot;       // permuted entry -> DAE entry
        std::vector<double> j_vals;      // permuted J values
        std::vector<double> c_vals;      // permuted C values
        bool prepared = false;
        int  n = 0;
        std::string error;
    } adj;

    // Builds the permuted pattern and the pivot-safe ordering.
    auto prepare = [&](const DaeMatrix& g, const DaeMatrix& c, double s) -> bool {
        const int n = static_cast<int>(g.rowPointers().size()) - 1;
        adj.n = n;
        const auto& rp = g.rowPointers();
        const auto& ci = g.colIndices();
        const auto& gv = g.values();
        const auto& cv = c.values();

        // The reference matrix for the ordering decision: J = G + s·C.
        std::vector<double> jref(gv.size(), 0.0);
        std::vector<double> zref(gv.size(), 0.0);
        for (size_t p = 0; p < gv.size(); ++p) {
            jref[p] = gv[p] + s * cv[p];
        }

        // Priority: branch rows make the best partners (their couplings are
        // the ±1 of a branch equation, so they give a well-scaled pivot).
        std::vector<int> priority(static_cast<size_t>(n), 0);
        for (int k = num_node_rows; k < n; ++k) {
            priority[static_cast<size_t>(k)] = 1;
        }

        StaticRowPermutation perm =
            StaticRowPermutation::build(rp, ci, jref, zref, priority);
        if (!perm.ok) {
            adj.error = "the transient Jacobian cannot be ordered for "
                        "factorization: " + perm.reason;
            return false;
        }

        adj.row_of_pos = perm.row_of_pos;
        adj.row_ptr.assign(static_cast<size_t>(n) + 1, 0);
        adj.col_idx.clear();
        adj.src_slot.clear();
        adj.j_vals.clear();
        adj.c_vals.clear();
        adj.col_idx.reserve(ci.size());
        adj.src_slot.reserve(ci.size());
        adj.j_vals.reserve(ci.size());
        adj.c_vals.reserve(ci.size());
        for (int k = 0; k < n; ++k) {
            adj.row_ptr[static_cast<size_t>(k)] = static_cast<int>(adj.col_idx.size());
            const int r = perm.row_of_pos[static_cast<size_t>(k)];
            for (int p = rp[r]; p < rp[r + 1]; ++p) {
                adj.col_idx.push_back(ci[static_cast<size_t>(p)]);
                adj.src_slot.push_back(p);
                adj.c_vals.push_back(cv[static_cast<size_t>(p)]);
                adj.j_vals.push_back(gv[static_cast<size_t>(p)] + s * cv[static_cast<size_t>(p)]);
            }
        }
        adj.row_ptr[static_cast<size_t>(n)] = static_cast<int>(adj.col_idx.size());

        adj.lu.setPattern(adj.row_ptr, adj.col_idx);
        adj.prepared = true;
        return true;
    };

    // Refills the values of the permuted matrix and factorizes it.
    auto refactorize = [&](const DaeMatrix& g, const DaeMatrix& c, double s) -> bool {
        const auto& gv = g.values();
        const auto& cv = c.values();
        for (size_t p = 0; p < adj.src_slot.size(); ++p) {
            const int src = adj.src_slot[p];
            adj.j_vals[p] = gv[static_cast<size_t>(src)] + s * cv[static_cast<size_t>(src)];
        }
        // A = G + s·C is REAL: it goes through the complex kernel with
        // omega = 0, so the imaginary half stays exactly zero.
        return adj.lu.factorize(/*omega=*/0.0, adj.j_vals.data(), adj.c_vals.data());
    };

    // ---- the sensitivity state ---------------------------------------------
    // Per parameter:
    //   y_prev[p]     the sensitivity vector at the previous accepted step
    //   dqdp_prev[p]  the explicit dependence of the devices' stored charges
    //                 at the previous point (the device hook reports it)
    //   e_prev[p]     the derivative of the stored dQ/dt (the trapezoidal
    //                 rule's secondary history)
    std::vector<std::vector<double>> y_prev(nparam);
    std::vector<std::vector<double>> dqdp_prev(nparam);
    std::vector<std::vector<double>> e_prev(nparam);
    std::vector<double> w_dae(static_cast<size_t>(n_dae), 0.0);
    std::vector<double> w_im(static_cast<size_t>(n_dae), 0.0);
    std::vector<double> dqdp_cur(static_cast<size_t>(n_dae), 0.0);
    std::vector<double> rhs(static_cast<size_t>(n_dae), 0.0);
    std::vector<double> rhs_im(static_cast<size_t>(n_dae), 0.0);
    std::vector<double> ct_y(static_cast<size_t>(n_dae), 0.0);   // C·y

    SensitivityContext ctx;
    ctx.kind = SensitivityContext::Kind::TRANSIENT;
    const bool trapezoidal = (method_ == IntegrationMethod::TRAPEZOIDAL);

    int failed_steps = 0;
    size_t step_index = 0;
    size_t kept_points = 0;
    static constexpr size_t kMaxRetainedPoints = 20000;

    // ---- the transient engine ----------------------------------------------
    TransientAnalysis tran(circuit_);
    tran.setIntegrationMethod(method_);
    tran.setMaxIterations(max_iterations_);
    tran.setTolerance(tolerance_);
    tran.setVerbose(verbose_);

    tran.setStepObserver([&](const TransientStep& info) -> bool {
        const DaeMatrix& g = *info.G;
        const DaeMatrix& c = *info.C;
        const int n = static_cast<int>(g.rowPointers().size()) - 1;

        // ---- factorize the effective matrix --------------------------------
        if (!adj.prepared) {
            if (!prepare(g, c, info.s)) return false;
        }
        const bool factored = refactorize(g, c, info.s);
        if (!factored) ++failed_steps;

        TransientSensitivityResult::Point pt;
        pt.time = info.time;
        pt.sens.assign(nparam, 0.0);
        pt.valid = factored;

        // ---- objective value ------------------------------------------------
        double obj = 0.0;
        const int n_use = std::min<int>(n, n_dae);
        for (int i = 0; i < n_use; ++i) {
            obj += c_dae[static_cast<size_t>(i)] *
                   (*info.x)[static_cast<size_t>(i)];
        }
        pt.obj_value = obj;

        if (factored) {
            const DaeSolution sol(*info.x, circuit_, node_map);
            ctx.s = info.s;
            ctx.time = info.time;

            for (size_t pi = 0; pi < nparam; ++pi) {
                if (!param_ok[pi]) continue;
                Device* dev = param_dev[pi];

                // The perturbation of the right-hand side, in the DAE space.
                std::fill(w_dae.begin(), w_dae.end(), 0.0);
                std::fill(w_im.begin(), w_im.end(), 0.0);
                std::fill(dqdp_cur.begin(), dqdp_cur.end(), 0.0);
                if (!dev->loadSensitivity(expanded[pi].second, sol, ctx,
                                          w_dae, w_im, dqdp_cur)) {
                    continue;
                }

                // HISTORY PROPAGATION: the accepted charge of the previous
                // step depends on the trajectory, hence on the parameter, and
                // that dependence is what carries the perturbation forward:
                //
                //   dQ_n/dp = (dQ/dx)·Y_n + dQ/dp|_explicit
                //
                // The trajectory part is the SAME for every device (every
                // device owns disjoint rows of dQ/dx), so the sum over the
                // devices is ONE sparse product with the charge Jacobian.
                const std::vector<double>& yp = y_prev[pi];
                const std::vector<double>& dqp = dqdp_prev[pi];
                if (!yp.empty() && ctx.s != 0.0) {
                    const auto& crp = c.rowPointers();
                    const auto& cci = c.colIndices();
                    const auto& ccv = c.values();
                    for (int i = 0; i < n_use; ++i) {
                        double sum = 0.0;
                        for (int p = crp[i]; p < crp[i + 1]; ++p) {
                            const int j = cci[static_cast<size_t>(p)];
                            if (j < 0 || j >= static_cast<int>(yp.size())) continue;
                            sum += ccv[static_cast<size_t>(p)] * yp[static_cast<size_t>(j)];
                        }
                        ct_y[static_cast<size_t>(i)] = sum;
                        w_dae[static_cast<size_t>(i)] +=
                            ctx.s * (sum + (i < static_cast<int>(dqp.size())
                                                ? dqp[static_cast<size_t>(i)] : 0.0));
                    }
                }

                // The trapezoidal rule's secondary history (the stored dQ/dt)
                // also depends on the parameter:
                //   e_{n+1} = s·(dQ_{n+1}/dp - dQ_n/dp) - e_n
                if (trapezoidal && !e_prev[pi].empty()) {
                    for (int i = 0; i < n_use; ++i) {
                        w_dae[static_cast<size_t>(i)] +=
                            e_prev[pi][static_cast<size_t>(i)];
                    }
                }

                // J·Y = w, solved on the PERMUTED system: the rows follow
                // the permutation.
                std::fill(rhs.begin(), rhs.end(), 0.0);
                std::fill(rhs_im.begin(), rhs_im.end(), 0.0);
                for (int k = 0; k < n_use; ++k) {
                    const int r = adj.row_of_pos[static_cast<size_t>(k)];
                    if (r < 0 || r >= n_dae) continue;
                    rhs[static_cast<size_t>(k)] = w_dae[static_cast<size_t>(r)];
                    rhs_im[static_cast<size_t>(k)] = w_im[static_cast<size_t>(r)];
                }
                adj.lu.solve(rhs, rhs_im);

                // Back to the DAE space, and the objective derivative.
                std::vector<double>& y = y_prev[pi];
                y.assign(static_cast<size_t>(n_use), 0.0);
                double sens = 0.0;
                for (int k = 0; k < n_use; ++k) {
                    const int r = adj.row_of_pos[static_cast<size_t>(k)];
                    if (r < 0 || r >= n_dae) continue;
                    y[static_cast<size_t>(r)] = rhs[static_cast<size_t>(k)];
                    sens += c_dae[static_cast<size_t>(r)] * rhs[static_cast<size_t>(k)];
                }
                pt.sens[pi] = sens;

                // ---- close the history recursion ----------------------------
                // dQ_{n+1}/dp = (dQ/dx)·Y_{n+1} + dQ/dp|_explicit
                // e_{n+1} = s·(dQ_{n+1}/dp - dQ_n/dp) - e_n
                if (trapezoidal) {
                    std::vector<double> e_new(static_cast<size_t>(n_use), 0.0);
                    const auto& crp = c.rowPointers();
                    const auto& cci = c.colIndices();
                    const auto& ccv = c.values();
                    for (int i = 0; i < n_use; ++i) {
                        double sum = 0.0;
                        for (int p = crp[i]; p < crp[i + 1]; ++p) {
                            const int j = cci[static_cast<size_t>(p)];
                            if (j < 0 || j >= static_cast<int>(y.size())) continue;
                            sum += ccv[static_cast<size_t>(p)] * y[static_cast<size_t>(j)];
                        }
                        const double dq_new = sum + dqdp_cur[static_cast<size_t>(i)];
                        const double dq_old = ct_y[static_cast<size_t>(i)] +
                            (i < static_cast<int>(dqp.size())
                                 ? dqp[static_cast<size_t>(i)] : 0.0);
                        const double e_old =
                            (i < static_cast<int>(e_prev[pi].size()))
                                ? e_prev[pi][static_cast<size_t>(i)] : 0.0;
                        e_new[static_cast<size_t>(i)] =
                            ctx.s * (dq_new - dq_old) - e_old;
                    }
                    e_prev[pi].swap(e_new);
                }
                dqdp_prev[pi] = dqdp_cur;
            }
        }

        if (csv_) {
            *csv_ << pt.time << "," << pt.obj_value;
            for (size_t pi = 0; pi < nparam; ++pi) *csv_ << "," << pt.sens[pi];
            *csv_ << "\n";
        }

        if (std::getenv("EXSPICE_DEBUG_SENS")) {
            std::fprintf(stderr, "[SENS-TRAN] t=%.6g obj=%.6g s=%.6g s0=%.6g s1=%.6g\n",
                         pt.time, pt.obj_value, info.s,
                         nparam > 0 ? pt.sens[0] : 0.0,
                         nparam > 1 ? pt.sens[1] : 0.0);
        }

        // The in-memory sample is decimated (the CSV keeps every point):
        // when the cap is reached the kept points are halved, so the
        // resolution degrades gracefully instead of the memory exploding.
        if (kept_points >= kMaxRetainedPoints) {
            std::vector<TransientSensitivityResult::Point> halved;
            halved.reserve(kept_points / 2 + 1);
            for (size_t i = 0; i < result.points.size(); i += 2) {
                halved.push_back(result.points[i]);
            }
            result.points.swap(halved);
            kept_points = result.points.size();
        }
        result.points.push_back(pt);
        ++kept_points;
        ++step_index;
        return true;
    });

    const TransientResult tran_result = tran.run(tstart_, tstop_, tstep_, use_mpi);
    if (!tran_result.success) {
        result.message = "transient simulation failed: " + tran_result.message;
        if (!adj.error.empty()) result.message += " (" + adj.error + ")";
        return result;
    }
    if (!adj.error.empty()) {
        result.message = adj.error;
        return result;
    }

    result.success = true;
    std::ostringstream msg;
    msg << step_index << " time step(s), " << n_dae
        << " unknowns, " << nparam << " parameter(s)";
    if (failed_steps > 0) {
        msg << ", " << failed_steps << " step(s) failed to factorize";
    }
    result.message = msg.str();
    return result;
}

} // namespace e_XSpice
