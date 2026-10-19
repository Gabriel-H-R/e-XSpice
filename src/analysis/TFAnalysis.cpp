/**
 * @author Gabriel Henrique Silva
 * @file TFAnalysis.cpp
 * @brief DC transfer function: OP + linearization + 2 solves.
 */

#include "TFAnalysis.h"

#include "../core/Device.h"
#include "../core/Node.h"
#include "../devices/Resistor.h"
#include "DaeSystem.h"
#include "../solver/ComplexSparseLU.h"
#include "SmallSignalLinearization.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

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

// Parse "V(node)" / "V(a,b)" / "I(dev)".
// Returns true on success; fills kind ('V' or 'I'), node1, node2 (for V(a,b)),
// device name (for I(dev)).
bool parseOutVar(const std::string& var, char& kind, std::string& a,
                 std::string& b, std::string& dev) {
    if (var.empty()) return false;
    std::string s = var;
    // strip surrounding whitespace
    size_t p = s.find_first_not_of(" \t");
    if (p == std::string::npos) return false;
    s = s.substr(p);
    p = s.find_last_not_of(" \t");
    if (p != std::string::npos) s = s.substr(0, p + 1);
    if (s.size() < 3 || s.back() != ')') return false;
    char k = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    if (k != 'V' && k != 'I') return false;
    if (s[1] != '(') return false;
    const std::string inner = s.substr(2, s.size() - 3);
    kind = k;
    if (k == 'V') {
        const size_t comma = inner.find(',');
        if (comma == std::string::npos) { a = inner; b.clear(); }
        else { a = inner.substr(0, comma); b = inner.substr(comma + 1); }
        dev.clear();
    } else {
        dev = inner;
        a.clear(); b.clear();
    }
    return true;
}

} // namespace

TFAnalysis::TFAnalysis(Circuit& circuit) : circuit_(circuit) {}
TFAnalysis::~TFAnalysis() = default;

TFResult TFAnalysis::run() {
    TFResult result;
    result.output_var = output_var_;
    result.input_source = input_source_;

    if (output_var_.empty()) {
        result.success = false;
        result.message = "no output variable given (use V(node) or I(dev))";
        return result;
    }
    if (input_source_.empty()) {
        result.success = false;
        result.message = "no input source given";
        return result;
    }

    char kind = 0;
    std::string a, b, dev;
    if (!parseOutVar(output_var_, kind, a, b, dev)) {
        result.success = false;
        result.message = "malformed output variable '" + output_var_ + "'";
        return result;
    }

    // ---- operating point + linearization (shared with .AC/.NOISE) ---------
    SmallSignalLinearization lin;
    {
        SmallSignalLinearization::Options opt;
        opt.max_iterations = max_iterations_;
        opt.abs_tol = std::max(1e-12, tolerance_ * 1e-3);
        opt.verbose = verbose_;
        opt.use_distributed_layout = false;     // .TF is single-rank
        if (!lin.prepare(circuit_, opt)) {
            result.success = false;
            result.message = lin.message();
            return result;
        }
    }
    const int nn = lin.size();
    const int n_full = lin.fullSize();

    // ---- resolve the input source ---------------------------------------
    Device* src = nullptr;
    for (const auto& d : circuit_.getDevices()) {
        if (iequals(d->getName(), input_source_)) { src = d.get(); break; }
    }
    if (!src) {
        result.success = false;
        result.message = "input source '" + input_source_ + "' not found";
        return result;
    }
    if (src->getType() != DeviceType::VOLTAGE_SOURCE &&
        src->getType() != DeviceType::CURRENT_SOURCE) {
        result.success = false;
        result.message = "input source '" + input_source_ +
                         "' is not an independent source";
        return result;
    }

    // ---- resolve the output variable ------------------------------------
    // Returns the frequency-system rows involved (pos, neg) and a kind tag:
    // 'V' node voltage, 'I' branch current of a device with a branch unknown.
    int out_pos = -1, out_neg = -1;
    bool out_is_current = false;
    if (kind == 'V') {
        const auto& node_map = lin.system().nodeMap();
        auto rowOfName = [&](const std::string& name, bool& found) -> int {
            found = false;
            if (name.empty() || name == "0") return -1;
            for (const auto& n : circuit_.getNodes()) {
                if (n->isGround()) continue;
                if (!iequals(n->getName(), name)) continue;
                const int id = n->getId();
                if (id < 0 || id >= static_cast<int>(node_map.size())) return -1;
                found = true;
                return lin.acRow(node_map[static_cast<size_t>(id)]);
            }
            return -1;
        };
        bool f1 = false, f2 = false;
        out_pos = rowOfName(a, f1);
        if (!f1) {
            result.success = false;
            result.message = "output node '" + a + "' not found";
            return result;
        }
        if (!b.empty()) {
            out_neg = rowOfName(b, f2);
            if (!f2) {
                result.success = false;
                result.message = "output node '" + b + "' not found";
                return result;
            }
        }
    } else {
        // I(dev): branch current of a voltage source / inductor, or the
        // conductance-weighted node difference for a resistor.
        Device* d = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), dev)) { d = dv.get(); break; }
        }
        if (!d) {
            result.success = false;
            result.message = "output device '" + dev + "' not found";
            return result;
        }
        if (d->getBranchIndex() >= 0) {
            // Voltage source / inductor: the branch unknown IS the current.
            out_pos = lin.acRow(d->getBranchIndex());
            out_neg = -2;     // sentinel: single-row value
        } else if (d->getType() == DeviceType::RESISTOR) {
            // Small-signal current = g · (V(pos) − V(neg)).
            out_pos = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(
                d->getNodes()[0]->getId())]);
            out_neg = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(
                d->getNodes()[1]->getId())]);
            out_is_current = true;   // handled specially below
        } else {
            result.success = false;
            result.message = "I(" + dev + ") is not a supported output "
                             "(use a voltage source, an inductor or a resistor)";
            return result;
        }
    }

    // ---- the matrix (already filtered + permuted) -------------------------
    const std::vector<double>& g = lin.gValues();
    const std::vector<double>& c = lin.cValues();
    const std::vector<int>& rp = lin.rowPointers();
    const std::vector<int>& ci = lin.colIndices();

    ComplexSparseLU lu;
    lu.setPattern(rp, ci);

    // ---- solve 1: drive the input source with a unit phasor ---------------
    std::vector<double> full_re(static_cast<size_t>(n_full), 0.0);
    std::vector<double> full_im(static_cast<size_t>(n_full), 0.0);
    if (src->getType() == DeviceType::VOLTAGE_SOURCE) {
        const int b = src->getBranchIndex();
        if (b >= 0) full_re[static_cast<size_t>(b)] += 1.0;
    } else {
        const auto& nodes = src->getNodes();
        const auto& nm = lin.system().nodeMap();
        if (nodes.size() >= 2 && nodes[0] && nodes[1]) {
            const int id0 = nodes[0]->getId();
            const int id1 = nodes[1]->getId();
            if (id0 >= 0 && id0 < static_cast<int>(nm.size())) {
                const int r = nm[static_cast<size_t>(id0)];
                if (r >= 0 && r < n_full) full_re[static_cast<size_t>(r)] -= 1.0;
            }
            if (id1 >= 0 && id1 < static_cast<int>(nm.size())) {
                const int r = nm[static_cast<size_t>(id1)];
                if (r >= 0 && r < n_full) full_re[static_cast<size_t>(r)] += 1.0;
            }
        }
    }
    std::vector<double> rhs_re = lin.permuteRhs(full_re);
    std::vector<double> rhs_im = lin.permuteRhs(full_im);
    std::vector<double> zero_im(static_cast<size_t>(nn), 0.0);

    if (std::getenv("EXSPICE_DEBUG_TF")) {
        std::fprintf(stderr, "[TF] nn=%d n_full=%d gnd=%d out_pos=%d out_neg=%d\n",
                     nn, n_full, lin.groundRow(), out_pos, out_neg);
        std::fprintf(stderr, "[TF] perm:");
        for (int k = 0; k < nn; ++k) {
            std::fprintf(stderr, " %d", lin.permutation().row_of_pos[static_cast<size_t>(k)]);
        }
        std::fprintf(stderr, "\n");
        std::fprintf(stderr, "[TF] matrix (G):\n");
        for (int i = 0; i < nn; ++i) {
            std::fprintf(stderr, "  row %d:", i);
            for (int p = rp[i]; p < rp[i + 1]; ++p) {
                std::fprintf(stderr, "  (%d,%.6g)", ci[p], g[p]);
            }
            std::fprintf(stderr, "\n");
        }
        // (rhs2/x2 are computed later; print them by running with a second pass)
    }
    if (!lu.factorize(0.0, g.data(), c.data())) {
        result.success = false;
        result.message = "could not factorize the DC system";
        return result;
    }
    std::vector<double> x_re = rhs_re, x_im = rhs_im;
    lu.solve(x_re, x_im);     // omega = 0 -> the factors are just G

    // Output variable from solve 1.
    double v_out = 0.0;
    if (kind == 'V') {
        if (out_pos >= 0) v_out += x_re[static_cast<size_t>(out_pos)];
        if (out_neg >= 0) v_out -= x_re[static_cast<size_t>(out_neg)];
    } else {
        // I(dev): branch current for V/L, or g·(V(pos)−V(neg)) for R.
        Device* d = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), dev)) { d = dv.get(); break; }
        }
        if (d && d->getBranchIndex() >= 0) {
            const int br = lin.acRow(d->getBranchIndex());
            if (br >= 0) v_out = x_re[static_cast<size_t>(br)];
        } else if (d && d->getType() == DeviceType::RESISTOR) {
            const double gdev = 1.0 / static_cast<const Resistor*>(d)->getResistance();
            const int p = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(
                d->getNodes()[0]->getId())]);
            const int n = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(
                d->getNodes()[1]->getId())]);
            const double vp = (p >= 0) ? x_re[static_cast<size_t>(p)] : 0.0;
            const double vn = (n >= 0) ? x_re[static_cast<size_t>(n)] : 0.0;
            v_out = gdev * (vp - vn);
        }
    }

    // Input resistance from solve 1.
    double rin = 0.0;
    if (src->getType() == DeviceType::VOLTAGE_SOURCE) {
        const int b = src->getBranchIndex();
        const int br = (b >= 0) ? lin.acRow(b) : -1;
        const double i_src = (br >= 0) ? x_re[static_cast<size_t>(br)] : 0.0;
        // i_src is the current ABSORBED by the source (the convention used by
        // the MNA stamp); the resistance the source sees is V/I_drawn, with
        // I_drawn = -i_src, hence the minus sign.
        rin = (std::fabs(i_src) > 1e-30) ? -1.0 / i_src : 0.0;
    } else {
        const auto& nodes = src->getNodes();
        int p = -1, n = -1;
        if (nodes.size() >= 2 && nodes[0] && nodes[1]) {
            const int id0 = nodes[0]->getId();
            const int id1 = nodes[1]->getId();
            if (id0 >= 0 && id0 < n_full) p = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(id0)]);
            if (id1 >= 0 && id1 < n_full) n = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(id1)]);
        }
        const double v_src = (p >= 0 ? x_re[static_cast<size_t>(p)] : 0.0)
                           - (n >= 0 ? x_re[static_cast<size_t>(n)] : 0.0);
        rin = v_src;     // I(src) = 1 A by construction
    }

    // ---- solve 2: inject a unit current at the output node ----------------
    // The input source stays in its small-signal state (short/open), which is
    // already what A encodes, so a single solve with a unit current at the
    // output node gives the output resistance directly.
    std::vector<double> full2(static_cast<size_t>(n_full), 0.0);
    if (kind == 'V') {
        if (out_pos >= 0) full2[static_cast<size_t>(lin.fullRow(out_pos))] += 1.0;
        if (out_neg >= 0) full2[static_cast<size_t>(lin.fullRow(out_neg))] -= 1.0;
    } else if (kind == 'I') {
        // For a current output, the "output node" for Rout is the device's
        // terminals — inject between them.
        Device* d = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), dev)) { d = dv.get(); break; }
        }
        if (d && d->getNodes().size() >= 2) {
            // Inject in the natural current direction of the device (pos -> neg)
            // so that Rout comes out positive for a passive load, matching the
            // V(out) convention where +1 goes into the output node.
            const auto& nm = lin.system().nodeMap();
            const int id0 = d->getNodes()[0]->getId();
            const int id1 = d->getNodes()[1]->getId();
            if (id0 >= 0 && id0 < static_cast<int>(nm.size())) {
                const int r = nm[static_cast<size_t>(id0)];
                if (r >= 0 && r < n_full) full2[static_cast<size_t>(r)] -= 1.0;
            }
            if (id1 >= 0 && id1 < static_cast<int>(nm.size())) {
                const int r = nm[static_cast<size_t>(id1)];
                if (r >= 0 && r < n_full) full2[static_cast<size_t>(r)] += 1.0;
            }
        }
    }
    std::vector<double> rhs2 = lin.permuteRhs(full2);
    std::vector<double> x2_re = rhs2, x2_im = zero_im;
    lu.solve(x2_re, x2_im);

    (void)rhs2;   // used only in the solve above

    double rout = 0.0;
    if (kind == 'V') {
        const double vp = (out_pos >= 0) ? x2_re[static_cast<size_t>(out_pos)] : 0.0;
        const double vn = (out_neg >= 0) ? x2_re[static_cast<size_t>(out_neg)] : 0.0;
        rout = vp - vn;
    } else {
        // For a current output, Rout is the resistance seen at the device's
        // terminals: inject 1 A and read V across them.
        Device* d = nullptr;
        for (const auto& dv : circuit_.getDevices()) {
            if (iequals(dv->getName(), dev)) { d = dv.get(); break; }
        }
        if (d && d->getNodes().size() >= 2) {
            const int id0 = d->getNodes()[0]->getId();
            const int id1 = d->getNodes()[1]->getId();
            const int p = (id0 >= 0 && id0 < n_full) ? lin.acRow(lin.system().nodeMap()[static_cast<size_t>(id0)]) : -1;
            const int n = (id1 >= 0 && id1 < n_full) ? lin.acRow(lin.system().nodeMap()[static_cast<size_t>(id1)]) : -1;
            const double vp = (p >= 0) ? x2_re[static_cast<size_t>(p)] : 0.0;
            const double vn = (n >= 0) ? x2_re[static_cast<size_t>(n)] : 0.0;
            // Report the magnitude: the sign of a resistance across two arbitrary
            // nodes is convention-dependent, and the magnitude is what a user
            // reads off a `.TF` line.
            rout = std::fabs(vp - vn);
        }
    }

    result.transfer_function = v_out;
    result.input_resistance  = rin;
    result.output_resistance  = rout;
    result.success = true;
    result.message = "transfer function = " + std::to_string(v_out) +
                     ", Rin = " + std::to_string(rin) +
                     " ohm, Rout = " + std::to_string(rout) + " ohm";
    return result;
}

} // namespace e_XSpice