/**
 * @author Gabriel Henrique Silva
 * @file ACSensitivityAnalysis.cpp
 * @brief AC sensitivity by the adjoint method.
 *
 * The analysis reuses the shared small-signal core (operating point, (G, C)
 * snapshot, pivot-safe ordering), so the AC sensitivity sees EXACTLY the
 * same system as the AC, noise, transfer-function and pole-zero analyses.
 *
 * Per frequency:
 *   1. factorize A(ω) = G + jωC (the pattern is frozen by the ordering);
 *   2. solve the AC system once to get the phasor X (the objective value);
 *   3. solve the ADJOINT system once to get λ = A^{-H}·(∂O/∂X);
 *   4. per parameter, ask the device for the right-hand-side perturbation w
 *      and take the scalar product λ^H·w.
 *
 * Step 4 costs a handful of flops per parameter, so the sweep cost is
 * independent of the number of parameters (the direct method would need one
 * extra solve per parameter per frequency).
 */

#include "ACSensitivityAnalysis.h"

#include "../core/DaeSolution.h"
#include "../core/Device.h"
#include "../core/Node.h"
#include "../parallel/MPIManager.h"
#include "../solver/ComplexSparseLU.h"
#include "SmallSignalLinearization.h"

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
constexpr double kPi2 = 6.28318530717958647692;   // 2*pi

bool iequals(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i]))) return false;
    }
    return true;
}
} // namespace

ACSensitivityAnalysis::ACSensitivityAnalysis(Circuit& circuit) : circuit_(circuit) {}
ACSensitivityAnalysis::~ACSensitivityAnalysis() = default;

ACSensitivityResult ACSensitivityAnalysis::run(bool use_mpi) {
    ACSensitivityResult result;
    result.variation = (variation_ == ACVariation::LIN)  ? "LIN"
                     : (variation_ == ACVariation::DEC)  ? "DEC" : "OCT";
    result.output_var = output_var_;

    const std::vector<double> grid =
        ACAnalysis::buildGrid(variation_, points_, fstart_, fstop_);
    if (grid.empty()) {
        result.message = "invalid frequency sweep specification";
        return result;
    }
    const size_t nfreq = grid.size();

    MPIManager& mpi = MPIManager::getInstance();
    const int rank = mpi.getRank();
    const int world = mpi.getSize();
    const bool distribute = point_distribution_ && use_mpi && world > 1 &&
                            nfreq >= static_cast<size_t>(2) * static_cast<size_t>(world);

    SmallSignalLinearization lin;
    {
        SmallSignalLinearization::Options opt;
        opt.max_iterations = max_iterations_;
        opt.abs_tol = std::max(1e-12, tolerance_ * 1e-3);
        opt.verbose = verbose_ && rank == 0;
        opt.use_distributed_layout = (use_mpi && !distribute);
        if (!lin.prepare(circuit_, opt)) {
            result.success = false;
            result.message = lin.message();
            return result;
        }
    }
    const int nn = lin.size();
    const int n_full = lin.fullSize();
    const auto& node_map = lin.system().nodeMap();

    // ---- output gradient (∂O/∂X) in DAE space ----------------------------
    std::vector<double> c_full(static_cast<size_t>(n_full), 0.0);
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
        const char kind = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
        const std::string inner = s.substr(2, s.size() - 3);

        if (kind == 'V') {
            const size_t comma = inner.find(',');
            const std::string na = (comma == std::string::npos) ? inner
                                                              : inner.substr(0, comma);
            const std::string nb = (comma == std::string::npos) ? std::string()
                                                              : inner.substr(comma + 1);
            auto rowOfName = [&](const std::string& name) -> int {
                if (name.empty() || name == "0") return -1;
                for (const auto& n : circuit_.getNodes()) {
                    if (n->isGround()) continue;
                    if (!iequals(n->getName(), name)) continue;
                    const int id = n->getId();
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
            if (r1 >= 0) c_full[static_cast<size_t>(r1)] += 1.0;
            if (r2 >= 0) c_full[static_cast<size_t>(r2)] -= 1.0;
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
                if (br >= 0 && br < n_full) c_full[static_cast<size_t>(br)] += 1.0;
            } else {
                result.message = "I(" + inner + ") is not a supported AC output";
                return result;
            }
        }
    }

    // The adjoint right-hand side lives in the ORIGINAL row indexing (the
    // transpose undoes the pivot-safe row permutation).
    const std::vector<double> c_orig = lin.dropGroundRow(c_full);
    // The objective value pairs the solution with the permuted selector.
    const std::vector<double> c_perm = lin.permuteRhs(c_full);

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

    // Resolve the devices once (the sweep must not search the device list
    // per parameter per frequency).
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
    for (size_t i = 0; i < nparam; ++i) {
        result.param_names.push_back(expanded[i].first + ":" + expanded[i].second);
    }

    // ---- per-thread LU slots ----------------------------------------------
    int nthreads = 1;
    #ifdef _OPENMP
    nthreads = std::max(1, omp_get_max_threads());
    #endif
    const int nslots = nthreads > 1 ? nthreads : 1;

    std::vector<ComplexSparseLU> lus(static_cast<size_t>(nslots));
    for (int s = 0; s < nslots; ++s) {
        lus[static_cast<size_t>(s)].setPattern(lin.rowPointers(), lin.colIndices());
    }

    // ---- per-thread scratch -----------------------------------------------
    // The frequencies are solved concurrently, so EVERY buffer is per slot:
    // the threads share nothing (no lock, no false sharing on the hot
    // vectors) and the result does not depend on the grouping.
    struct Scratch {
        std::vector<double> x_re, x_im;
        std::vector<double> lam_re, lam_im;
        std::vector<double> ac_re, ac_im;
        std::vector<double> w_re, w_im;
        std::vector<double> dqdp;
        std::vector<double> x_dae_re, x_dae_im;
    };
    std::vector<Scratch> scratch(static_cast<size_t>(nslots));
    for (auto& s : scratch) {
        s.x_re.assign(static_cast<size_t>(nn), 0.0);
        s.x_im.assign(static_cast<size_t>(nn), 0.0);
        s.lam_re.assign(static_cast<size_t>(nn), 0.0);
        s.lam_im.assign(static_cast<size_t>(nn), 0.0);
        s.ac_re.assign(static_cast<size_t>(n_full), 0.0);
        s.ac_im.assign(static_cast<size_t>(n_full), 0.0);
        s.w_re.assign(static_cast<size_t>(n_full), 0.0);
        s.w_im.assign(static_cast<size_t>(n_full), 0.0);
        s.dqdp.assign(static_cast<size_t>(n_full), 0.0);
        s.x_dae_re.assign(static_cast<size_t>(n_full), 0.0);
        s.x_dae_im.assign(static_cast<size_t>(n_full), 0.0);
    }

    // ---- CSV header --------------------------------------------------------
    if (csv_) {
        csv_->precision(17);
        if (!csv_header_written_) {
            *csv_ << "frequency,obj_mag,obj_phase_deg";
            for (const auto& name : result.param_names) {
                *csv_ << "," << name << "_mag," << name << "_phase_deg";
            }
            *csv_ << "\n";
            csv_header_written_ = true;
        }
    }

    // ---- retention ---------------------------------------------------------
    size_t stride = 1;
    if (max_retained_ > 0 && nfreq > max_retained_) {
        stride = (nfreq + max_retained_ - 1) / max_retained_;
    }
    const auto retain = [this, stride, nfreq](size_t index) -> bool {
        return max_retained_ == 0 || (index % stride == 0) || (index + 1 == nfreq);
    };

    long long flops_total = 0;
    int failures_total = 0;
    (void)flops_total;   // the factorization cost is per-thread, not additive

    // ---- the sweep ---------------------------------------------------------
    // Frequencies are INDEPENDENT, so they are processed concurrently on the
    // replicated (G, C): each thread owns one ComplexSparseLU (its own factor
    // arrays) plus its own buffers — no sharing, no lock, and the result does
    // not depend on how the frequencies were grouped. The sweep runs in BLOCKS
    // so the RAM stays bounded on a long grid.
    const size_t payload_cap = static_cast<size_t>(1) << 20;
    const size_t block = std::max<size_t>(
        1, std::min<size_t>(4096,
                             payload_cap / std::max<size_t>(
                                 2, static_cast<size_t>(nn) *
                                        static_cast<size_t>(nparam + 4))));

    for (size_t base = 0; base < nfreq; base += block) {
        const size_t hi = std::min(nfreq, base + block);
        const int count = static_cast<int>(hi - base);

        std::vector<ACSensitivityResult::Point> pts(static_cast<size_t>(count));

        #pragma omp parallel for schedule(dynamic, 1)
        for (int k = 0; k < count; ++k) {
            const size_t fi = base + static_cast<size_t>(k);
            const double omega = kPi2 * grid[fi];
            #ifdef _OPENMP
            const int slot = omp_get_thread_num() % nslots;
            #else
            const int slot = 0;
            #endif
            ComplexSparseLU& lu = lus[static_cast<size_t>(slot)];
            Scratch& sc = scratch[static_cast<size_t>(slot)];
            std::vector<double>& x_re = sc.x_re;
            std::vector<double>& x_im = sc.x_im;
            std::vector<double>& lam_re = sc.lam_re;
            std::vector<double>& lam_im = sc.lam_im;
            std::vector<double>& ac_re = sc.ac_re;
            std::vector<double>& ac_im = sc.ac_im;
            std::vector<double>& w_re = sc.w_re;
            std::vector<double>& w_im = sc.w_im;
            std::vector<double>& dqdp_scratch = sc.dqdp;
            ACSensitivityResult::Point& pt = pts[static_cast<size_t>(k)];
            pt.frequency = grid[fi];
            pt.sens_mag.assign(nparam, 0.0);
            pt.sens_phase_deg.assign(nparam, 0.0);

            // The context is per iteration (the threads share nothing).
            SensitivityContext ctx;
            ctx.kind = SensitivityContext::Kind::AC;
            ctx.omega = omega;

            if (!lu.factorize(omega, lin.gValues().data(), lin.cValues().data())) {
                pt.valid = false;
                continue;
            }

            // ---- AC excitation and solution X ------------------------------
            std::fill(ac_re.begin(), ac_re.end(), 0.0);
            std::fill(ac_im.begin(), ac_im.end(), 0.0);
            for (const auto& d : circuit_.getDevices()) {
                d->loadACExcitation(omega, ac_re, ac_im);
            }
            // Rows follow the permutation; columns stay in the DAE indexing.
            for (int i = 0; i < nn; ++i) {
                const int r = lin.fullRow(i);
                if (r >= 0 && r < n_full) {
                    x_re[static_cast<size_t>(i)] = ac_re[static_cast<size_t>(r)];
                    x_im[static_cast<size_t>(i)] = ac_im[static_cast<size_t>(r)];
                }
            }
            lu.solve(x_re, x_im);

            // ---- adjoint: A^H·λ = c ----------------------------------------
            // A^H·λ = c  <=>  A^T·conj(λ) = conj(c): the transposed solve
            // takes conj(c) in the ORIGINAL row indexing and returns conj(λ)
            // in the permuted indexing.
            for (int i = 0; i < nn; ++i) {
                lam_re[static_cast<size_t>(i)] = c_orig[static_cast<size_t>(i)];
                lam_im[static_cast<size_t>(i)] = 0.0;
            }
            lu.solveTranspose(lam_re, lam_im);
            for (int i = 0; i < nn; ++i) {
                lam_im[static_cast<size_t>(i)] = -lam_im[static_cast<size_t>(i)];
            }

            // ---- objective value --------------------------------------------
            {
                double re = 0.0, im = 0.0;
                for (int i = 0; i < nn; ++i) {
                    re += c_perm[static_cast<size_t>(i)] * x_re[static_cast<size_t>(i)];
                    im += c_perm[static_cast<size_t>(i)] * x_im[static_cast<size_t>(i)];
                }
                pt.obj_mag = std::hypot(re, im);
                pt.obj_phase_deg =
                    std::atan2(im, re) * 180.0 / 3.14159265358979323846;
            }

            // ---- per-parameter sensitivity ----------------------------------
            // The device reads the PHASOR solution, so it needs the DAE-space
            // complex vector.
            lin.unpermuteSolution(x_re, sc.x_dae_re);
            lin.unpermuteSolution(x_im, sc.x_dae_im);
            const DaeSolution sol(sc.x_dae_re, sc.x_dae_im, circuit_, node_map);

            for (size_t pi = 0; pi < nparam; ++pi) {
                if (!param_ok[pi]) continue;
                Device* dev = param_dev[pi];
                std::fill(w_re.begin(), w_re.end(), 0.0);
                std::fill(w_im.begin(), w_im.end(), 0.0);
                if (!dev->loadSensitivity(expanded[pi].second, sol, ctx,
                                          w_re, w_im, dqdp_scratch)) {
                    continue;
                }
                const std::vector<double> wp_re = lin.permuteRhs(w_re);
                const std::vector<double> wp_im = lin.permuteRhs(w_im);

                // dO/dp = λ^H·w
                double re = 0.0, im = 0.0;
                for (int i = 0; i < nn; ++i) {
                    const double lr = lam_re[static_cast<size_t>(i)];
                    const double li = lam_im[static_cast<size_t>(i)];
                    const double vr = wp_re[static_cast<size_t>(i)];
                    const double vi = wp_im[static_cast<size_t>(i)];
                    // conj(λ)·w = (lr - j·li)(vr + j·vi)
                    re += lr * vr + li * vi;
                    im += lr * vi - li * vr;
                }
                pt.sens_mag[pi] = std::hypot(re, im);
                pt.sens_phase_deg[pi] =
                    std::atan2(im, re) * 180.0 / 3.14159265358979323846;
            }
        }

        // ---- collect the block (in order) ----------------------------------
        for (int k = 0; k < count; ++k) {
            const size_t fi = base + static_cast<size_t>(k);
            ACSensitivityResult::Point& pt = pts[static_cast<size_t>(k)];
            if (!pt.valid) ++failures_total;

            if (csv_) {
                *csv_ << pt.frequency << "," << pt.obj_mag << ","
                      << pt.obj_phase_deg;
                for (size_t pi = 0; pi < nparam; ++pi) {
                    *csv_ << "," << pt.sens_mag[pi] << ","
                          << pt.sens_phase_deg[pi];
                }
                *csv_ << "\n";
            }
            if (retain(fi)) {
                result.points.push_back(std::move(pt));
            }
        }
    }

    result.success = true;
    std::ostringstream msg;
    msg << nfreq << " frequencies, " << nn << " unknowns, "
        << nparam << " parameter(s)";
    if (failures_total > 0) msg << ", " << failures_total << " frequency point(s) failed";
    result.message = msg.str();
    return result;
}

} // namespace e_XSpice
