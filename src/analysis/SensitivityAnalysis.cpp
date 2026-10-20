/**
 * @author Gabriel Henrique Silva
 * @file SensitivityAnalysis.cpp
 * @brief DC sensitivity by the adjoint method.
 *
 * The analysis reuses the shared small-signal core (operating point, (G, C)
 * snapshot, pivot-safe ordering) so that the DC sensitivity sees EXACTLY the
 * same system as the AC, noise, transfer-function and pole-zero analyses.
 * The adjoint solve is SPARSE (one LU factorization of G + a transposed
 * triangular sweep), not the dense O(n^3) solve an earlier version used.
 */

#include "SensitivityAnalysis.h"

#include "../core/DaeSolution.h"
#include "../core/Device.h"
#include "../core/Node.h"
#include "../devices/Resistor.h"
#include "../solver/ComplexSparseLU.h"
#include "DaeSystem.h"
#include "SmallSignalLinearization.h"

#include <algorithm>
#include <cmath>
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

SensitivityAnalysis::SensitivityAnalysis(Circuit& circuit) : circuit_(circuit) {}
SensitivityAnalysis::~SensitivityAnalysis() = default;

SensitivityResult SensitivityAnalysis::run() {
    SensitivityResult result;
    result.output_var = output_var_;

    if (output_var_.empty()) {
        result.success = false;
        result.message = "no output variable given (use V(node) or I(dev))";
        return result;
    }

    // ---- operating point + linearization ----------------------------------
    // The shared core runs the OP, snapshots (G, C), drops the ground row
    // and freezes a pivot-safe ordering. reference_omega = 0 because the DC
    // matrix is G ALONE: an inductor's branch column has no diagonal there,
    // and pretending otherwise (j·ω·L) would pick a pivot that does not
    // exist.
    SmallSignalLinearization lin;
    {
        SmallSignalLinearization::Options opt;
        opt.max_iterations = 200;
        opt.abs_tol = 1e-12;
        opt.verbose = verbose_;
        opt.reference_omega = 0.0;
        if (!lin.prepare(circuit_, opt)) {
            result.success = false;
            result.message = lin.message();
            return result;
        }
    }

    const int nn = lin.size();
    const int n_full = lin.fullSize();
    const auto& node_map = lin.system().nodeMap();

    // ---- output gradient (∂O/∂x) in DAE space ----------------------------
    std::vector<double> c_full(static_cast<size_t>(n_full), 0.0);
    {
        std::string s = output_var_;
        size_t p = s.find_first_not_of(" \t");
        if (p != std::string::npos) s = s.substr(p);
        p = s.find_last_not_of(" \t");
        if (p != std::string::npos) s = s.substr(0, p + 1);
        if (s.size() < 3 || s.back() != ')' || s[1] != '(') {
            result.success = false;
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
                result.success = false;
                result.message = "output node of '" + output_var_ + "' not found";
                return result;
            }
            if (r1 >= 0) c_full[static_cast<size_t>(r1)] += 1.0;
            if (r2 >= 0) c_full[static_cast<size_t>(r2)] -= 1.0;
        } else {
            // I(dev): the branch current of a voltage-defined device, or the
            // terminal current of a resistor.
            Device* d = nullptr;
            for (const auto& dv : circuit_.getDevices()) {
                if (iequals(dv->getName(), inner)) { d = dv.get(); break; }
            }
            if (!d) {
                result.success = false;
                result.message = "output device '" + inner + "' not found";
                return result;
            }
            if (d->getBranchIndex() >= 0) {
                const int br = d->getBranchIndex();
                if (br >= 0 && br < n_full) c_full[static_cast<size_t>(br)] += 1.0;
            } else if (d->getType() == DeviceType::RESISTOR) {
                // I = (v+ - v-)/R  =>  dI/dv+ = 1/R, dI/dv- = -1/R.
                const double gdev = 1.0 / static_cast<const Resistor*>(d)->getResistance();
                const Node* np = d->getNodes()[0];
                const Node* nn = d->getNodes()[1];
                if (np && !np->isGround()) {
                    const int id = np->getId();
                    if (id >= 0 && id < static_cast<int>(node_map.size())) {
                        const int r = node_map[static_cast<size_t>(id)];
                        if (r >= 0) c_full[static_cast<size_t>(r)] += gdev;
                    }
                }
                if (nn && !nn->isGround()) {
                    const int id = nn->getId();
                    if (id >= 0 && id < static_cast<int>(node_map.size())) {
                        const int r = node_map[static_cast<size_t>(id)];
                        if (r >= 0) c_full[static_cast<size_t>(r)] -= gdev;
                    }
                }
            } else {
                result.success = false;
                result.message = "I(" + inner + ") is not a supported output";
                return result;
            }
        }
    }

    // ---- adjoint solve: A^T·λ = c ----------------------------------------
    // The TRANSPOSED solve takes its right-hand side in the ORIGINAL row
    // indexing (the transpose undoes the pivot-safe row permutation) and
    // returns the solution in the permuted indexing.
    ComplexSparseLU lu;
    lu.setPattern(lin.rowPointers(), lin.colIndices());
    if (!lu.factorize(/*omega=*/0.0, lin.gValues().data(), lin.cValues().data())) {
        result.success = false;
        result.message = "the operating-point Jacobian could not be factorized: " +
                         lu.message();
        return result;
    }
    if (lu.perturbedPivot() && verbose_) {
        std::cout << "[SENS] warning: a pivot had to be perturbed "
                     "(the operating point is singular for this parameter "
                     "combination)" << std::endl;
    }

    std::vector<double> lam_re = lin.dropGroundRow(c_full);
    std::vector<double> lam_im(static_cast<size_t>(nn), 0.0);
    lu.solveTranspose(lam_re, lam_im);

    // ---- parameter list ----------------------------------------------------
    // An empty parameter name means "every sensitivity parameter of the
    // device" (that is what `.SENS V(out) R1` means).
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

    // ---- per-parameter sensitivity ----------------------------------------
    SensitivityContext ctx;
    ctx.kind = SensitivityContext::Kind::DC;
    ctx.time = 0.0;

    const DaeSolution sol(lin.system().solution(), circuit_, node_map);
    std::vector<double> w_re(static_cast<size_t>(n_full), 0.0);
    std::vector<double> w_im(static_cast<size_t>(n_full), 0.0);
    std::vector<double> dqdp_scratch(static_cast<size_t>(n_full), 0.0);

    for (const auto& [dev_name, param_name] : expanded) {
        result.param_names.push_back(dev_name + ":" + param_name);

        Device* dev = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), dev_name)) { dev = dv.get(); break; }
        }
        if (!dev) {
            result.sensitivities.push_back(0.0);
            continue;
        }
        bool supported = false;
        for (const auto& pn : dev->sensitivityParameters()) {
            if (iequals(pn, param_name)) { supported = true; break; }
        }
        if (!supported) {
            result.sensitivities.push_back(0.0);
            continue;
        }

        std::fill(w_re.begin(), w_re.end(), 0.0);
        std::fill(w_im.begin(), w_im.end(), 0.0);
        if (!dev->loadSensitivity(param_name, sol, ctx, w_re, w_im, dqdp_scratch)) {
            result.sensitivities.push_back(0.0);
            continue;
        }

        // λ lives in the permuted indexing; w has to follow the rows.
        const std::vector<double> w_perm = lin.permuteRhs(w_re);
        double sens = 0.0;
        for (int i = 0; i < nn; ++i) {
            sens += lam_re[static_cast<size_t>(i)] * w_perm[static_cast<size_t>(i)];
        }
        result.sensitivities.push_back(sens);
    }

    result.success = true;
    std::ostringstream msg;
    msg << result.param_names.size() << " parameter(s), "
        << "output: " << output_var_
        << ", " << nn << " unknowns";
    result.message = msg.str();

    if (csv_) {
        csv_->precision(17);
        *csv_ << "parameter,sensitivity\n";
        for (size_t i = 0; i < result.param_names.size(); ++i) {
            *csv_ << result.param_names[i] << "," << result.sensitivities[i] << "\n";
        }
    }
    return result;
}

} // namespace e_XSpice
