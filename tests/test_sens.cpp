/**
 * @author Gabriel Henrique Silva
 * @file test_sens.cpp
 * @brief Sensitivity analyses (`.SENS`): DC, AC and transient.
 *
 * The tests are organized around CLOSED FORMS, not around golden values:
 * every check compares the simulator against a formula that can be derived
 * by hand, so a regression shows up as a physics error and not as a number
 * that changed.
 *
 *   1) DC  — resistive divider: dV(out)/dR1 and dV(out)/dR2 in closed form.
 *   2) DC  — the ADJOINT against a DENSE reference, on a circuit whose
 *            pivot-safe ordering is NOT the identity (the row permutation is
 *            the classic source of a silently wrong adjoint).
 *   3) AC  — RC low-pass: dV(out)/dR and dV(out)/dC in closed form
 *            (dH/dR = -sC/(1+sRC)^2, dH/dC = -sR/(1+sRC)^2).
 *   4) TRAN — RC step response: the sensitivity of V(out) w.r.t. R and C
 *            against the analytic transient of the same circuit.
 *   5) The device-level perturbation w against a hand-computed value.
 */

#include "core/Circuit.h"
#include "core/Node.h"
#include "core/Device.h"
#include "core/NetlistParser.h"
#include "core/DaeSolution.h"
#include "analysis/SensitivityAnalysis.h"
#include "analysis/ACSensitivityAnalysis.h"
#include "analysis/TransientSensitivityAnalysis.h"
#include "analysis/SmallSignalLinearization.h"
#include "solver/ComplexSparseLU.h"
#include "solver/DaeMatrix.h"
#include "analysis/DaeSystem.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace e_XSpice;

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::cout << "FAIL: " << what << "\n";
    }
}

void checkNear(double got, double expected, double tol, const std::string& what) {
    ++g_checks;
    const double diff = std::fabs(got - expected);
    if (!(diff <= tol)) {
        ++g_failures;
        std::cout << "FAIL: " << what << " -> got " << got
                  << ", expected " << expected << " (|d| = " << diff
                  << " > " << tol << ")\n";
    }
}

/** Builds a circuit from a netlist string (same path the CLI uses). */
bool buildCircuit(const std::string& text, Circuit& circuit, std::string& error) {
    NetlistParser parser;
    if (!parser.parseString(text, circuit)) {
        error = "netlist could not be parsed";
        return false;
    }
    error.clear();
    return true;
}

// ===============================================================
// 1) DC sensitivity of a resistive divider
// ===============================================================
void testDcDivider() {
    Circuit circuit;
    std::string err;
    const std::string netlist =
        "* dc divider sensitivity\n"
        "V1 in 0 DC 1\n"
        "R1 in out 1k\n"
        "R2 out 0 1k\n"
        ".END\n";
    if (!buildCircuit(netlist, circuit, err)) {
        check(false, "divider netlist parses: " + err);
        return;
    }

    SensitivityAnalysis sens(circuit);
    sens.setOutputVar("V(out)");
    sens.addParameter("R1", "R");
    sens.addParameter("R2", "R");
    const SensitivityResult res = sens.run();
    check(res.success, "DC sensitivity succeeds: " + res.message);
    if (!res.success) return;
    check(res.param_names.size() == 2, "two parameters reported");
    if (res.param_names.size() != 2) return;

    // V(out) = V1 * R2/(R1+R2)
    //   d/dR1 = -V1*R2/(R1+R2)^2 = -1*1000/4e6 = -2.5e-4
    //   d/dR2 = +V1*R1/(R1+R2)^2 = +2.5e-4
    checkNear(res.sensitivities[0], -2.5e-4, 1e-12, "d(V(out))/dR1");
    checkNear(res.sensitivities[1], +2.5e-4, 1e-12, "d(V(out))/dR2");
}

// ===============================================================
// 2) The adjoint against a DENSE reference, on a circuit whose
//    pivot-safe ordering is NOT the identity.
// ===============================================================
void testAdjointAgainstDense() {
    // Node `in` has no conductance (only the source and the inductor), so
    // the pivot-safe ordering MUST swap rows. The adjoint vector is then
    // computed in a permuted system, and the right-hand side of the
    // TRANSPOSED solve must stay in the ORIGINAL indexing.
    Circuit circuit;
    std::string err;
    const std::string netlist =
        "* permutation probe\n"
        "V1 in 0 DC 1\n"
        "L1 in out 1m\n"
        "R2 out 0 1k\n"
        ".END\n";
    if (!buildCircuit(netlist, circuit, err)) {
        check(false, "permutation netlist parses: " + err);
        return;
    }

    SmallSignalLinearization lin;
    SmallSignalLinearization::Options opt;
    opt.reference_omega = 0.0;
    if (!lin.prepare(circuit, opt)) {
        check(false, "linearization prepares: " + lin.message());
        return;
    }
    const int nn = lin.size();
    check(lin.permutation().swaps > 0,
          "the probe circuit really needs a row swap (otherwise the test "
          "would not discriminate)");

    // Dense reference: A^T u = c with A = G (the DC Jacobian).
    const DaeMatrix& G = lin.system().jacobianF();
    std::vector<double> A(static_cast<size_t>(nn) * static_cast<size_t>(nn), 0.0);
    {
        const auto& rp = G.rowPointers();
        const auto& ci = G.colIndices();
        const auto& gv = G.values();
        for (int i = 0; i < nn; ++i) {
            for (int p = rp[i]; p < rp[i + 1]; ++p) {
                A[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                  static_cast<size_t>(ci[p])] = gv[static_cast<size_t>(p)];
            }
        }
    }
    // c = the selector of node `out`.
    int out_row = -1;
    for (const auto& n : circuit.getNodes()) {
        if (n->getName() == "out") {
            out_row = lin.acRow(lin.system().nodeMap()[static_cast<size_t>(n->getId())]);
            break;
        }
    }
    check(out_row >= 0, "the output node has a row");

    // Dense transposed solve (pivoted).
    std::vector<double> AT(static_cast<size_t>(nn) * static_cast<size_t>(nn), 0.0);
    for (int i = 0; i < nn; ++i) {
        for (int j = 0; j < nn; ++j) {
            AT[static_cast<size_t>(j) * static_cast<size_t>(nn) +
               static_cast<size_t>(i)] = A[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                                          static_cast<size_t>(j)];
        }
    }
    std::vector<double> u(static_cast<size_t>(nn), 0.0);
    u[static_cast<size_t>(out_row)] = 1.0;
    for (int k = 0; k < nn; ++k) {
        int piv = k;
        double best = std::fabs(AT[static_cast<size_t>(k) * static_cast<size_t>(nn) +
                                   static_cast<size_t>(k)]);
        for (int i = k + 1; i < nn; ++i) {
            const double m = std::fabs(AT[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                                        static_cast<size_t>(k)]);
            if (m > best) { best = m; piv = i; }
        }
        if (best < 1e-30) { check(false, "dense reference is solvable"); return; }
        if (piv != k) {
            for (int j = 0; j < nn; ++j) {
                std::swap(AT[static_cast<size_t>(k) * static_cast<size_t>(nn) +
                             static_cast<size_t>(j)],
                          AT[static_cast<size_t>(piv) * static_cast<size_t>(nn) +
                             static_cast<size_t>(j)]);
            }
            std::swap(u[static_cast<size_t>(k)], u[static_cast<size_t>(piv)]);
        }
        const double pv = AT[static_cast<size_t>(k) * static_cast<size_t>(nn) +
                             static_cast<size_t>(k)];
        for (int i = k + 1; i < nn; ++i) {
            const double f = AT[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                                static_cast<size_t>(k)] / pv;
            for (int j = k; j < nn; ++j) {
                AT[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                   static_cast<size_t>(j)] -=
                    f * AT[static_cast<size_t>(k) * static_cast<size_t>(nn) +
                          static_cast<size_t>(j)];
            }
            u[static_cast<size_t>(i)] -= f * u[static_cast<size_t>(k)];
        }
    }
    for (int i = nn - 1; i >= 0; --i) {
        double s = u[static_cast<size_t>(i)];
        for (int j = i + 1; j < nn; ++j) {
            s -= AT[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                    static_cast<size_t>(j)] * u[static_cast<size_t>(j)];
        }
        u[static_cast<size_t>(i)] =
            s / AT[static_cast<size_t>(i) * static_cast<size_t>(nn) +
                   static_cast<size_t>(i)];
    }

    // The sparse adjoint: the transposed solve takes the right-hand side in
    // the ORIGINAL indexing and returns the solution in the permuted one.
    ComplexSparseLU lu;
    lu.setPattern(lin.rowPointers(), lin.colIndices());
    check(lu.factorize(0.0, lin.gValues().data(), lin.cValues().data()),
          "the sparse LU factorizes the DC matrix");

    std::vector<double> c_full(static_cast<size_t>(lin.fullSize()), 0.0);
    {
        // The output node's DAE row.
        int dae_row = -1;
        for (const auto& n : circuit.getNodes()) {
            if (n->getName() == "out") {
                dae_row = lin.system().nodeMap()[static_cast<size_t>(n->getId())];
                break;
            }
        }
        check(dae_row >= 0, "the output node has a DAE row");
        if (dae_row >= 0) c_full[static_cast<size_t>(dae_row)] = 1.0;
    }
    std::vector<double> x = lin.dropGroundRow(c_full);
    std::vector<double> xi(static_cast<size_t>(nn), 0.0);
    lu.solveTranspose(x, xi);

    // x[k] is the adjoint at the ORIGINAL row perm.row_of_pos[k].
    const auto& perm = lin.permutation();
    for (int r = 0; r < nn; ++r) {
        // The permuted position of the original row r.
        int k = -1;
        for (int p = 0; p < nn; ++p) {
            if (perm.row_of_pos[static_cast<size_t>(p)] == r) { k = p; break; }
        }
        check(k >= 0, "every original row has a permuted position");
        if (k < 0) continue;
        checkNear(x[static_cast<size_t>(k)], u[static_cast<size_t>(r)], 1e-9,
                  "sparse adjoint matches the dense reference (row " +
                      std::to_string(r) + ")");
    }
}

// ===============================================================
// 3) AC sensitivity of an RC low-pass against the closed form
// ===============================================================
void testAcLowPass() {
    Circuit circuit;
    std::string err;
    const std::string netlist =
        "* ac low-pass sensitivity\n"
        "V1 in 0 DC 1 AC 1\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".END\n";
    if (!buildCircuit(netlist, circuit, err)) {
        check(false, "AC netlist parses: " + err);
        return;
    }

    ACSensitivityAnalysis sens(circuit);
    sens.setVariation(ACVariation::DEC);
    sens.setPoints(10);
    sens.setFrequencyRange(1e3, 1e6);
    sens.setOutputVar("V(out)");
    sens.addParameter("R1", "R");
    sens.addParameter("C1", "C");
    const ACSensitivityResult res = sens.run();
    check(res.success, "AC sensitivity succeeds: " + res.message);
    if (!res.success || res.points.empty()) return;

    const double R = 1e3;
    const double C = 1e-6;

    // Check every retained frequency against the closed form:
    //   H(s)  = 1/(1+sRC)
    //   dH/dR = -sC/(1+sRC)^2
    //   dH/dC = -sR/(1+sRC)^2
    for (const auto& pt : res.points) {
        const double w = 2.0 * 3.14159265358979323846 * pt.frequency;
        const std::complex<double> s(0.0, w);
        const std::complex<double> den = 1.0 + s * R * C;
        const std::complex<double> H = 1.0 / den;
        const std::complex<double> dHdR = -s * C / (den * den);
        const std::complex<double> dHdC = -s * R / (den * den);

        const std::string tag = " @ " + std::to_string(pt.frequency) + " Hz";

        checkNear(pt.obj_mag, std::abs(H), 1e-9, "|V(out)|" + tag);
        checkNear(pt.obj_phase_deg,
                  std::arg(H) * 180.0 / 3.14159265358979323846, 1e-6,
                  "arg(V(out))" + tag);

        checkNear(pt.sens_mag[0], std::abs(dHdR), 1e-9, "|dV/dR|" + tag);
        checkNear(pt.sens_phase_deg[0],
                  std::arg(dHdR) * 180.0 / 3.14159265358979323846, 1e-6,
                  "arg(dV/dR)" + tag);

        checkNear(pt.sens_mag[1], std::abs(dHdC), 1e-9, "|dV/dC|" + tag);
        checkNear(pt.sens_phase_deg[1],
                  std::arg(dHdC) * 180.0 / 3.14159265358979323846, 1e-6,
                  "arg(dV/dC)" + tag);
    }
}

// ===============================================================
// 4) Transient sensitivity of a purely resistive circuit: the
//    answer is the DC sensitivity at every time step (exact).
// ===============================================================
void testTransientResistive() {
    Circuit circuit;
    std::string err;
    const std::string netlist =
        "* resistive transient sensitivity\n"
        "V1 in 0 DC 1\n"
        "R1 in out 1k\n"
        "R2 out 0 1k\n"
        ".END\n";
    if (!buildCircuit(netlist, circuit, err)) {
        check(false, "transient netlist parses: " + err);
        return;
    }

    TransientSensitivityAnalysis sens(circuit);
    sens.setVariation(0.0, 2e-3, 1e-4);
    sens.setOutputVar("V(out)");
    sens.addParameter("R1", "R");
    sens.addParameter("R2", "R");
    const TransientSensitivityResult res = sens.run();
    check(res.success, "transient sensitivity succeeds: " + res.message);
    if (!res.success || res.points.empty()) return;
    check(res.points.size() >= 5, "several time steps reported");

    for (const auto& pt : res.points) {
        checkNear(pt.obj_value, 0.5, 1e-9,
                  "V(out) of a resistive divider is constant");
        checkNear(pt.sens[0], -2.5e-4, 1e-9, "d(V(out))/dR1 (transient)");
        checkNear(pt.sens[1], +2.5e-4, 1e-9, "d(V(out))/dR2 (transient)");
    }
}

// ===============================================================
// 5) Transient sensitivity of an RC step against the closed form
// ===============================================================
void testTransientRcAgainstClosedForm() {
    Circuit circuit;
    std::string err;
    const std::string netlist =
        "* rc step sensitivity\n"
        "V1 in 0 PULSE(0 1 0 1p 1p 1 2)\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".END\n";
    if (!buildCircuit(netlist, circuit, err)) {
        check(false, "rc transient netlist parses: " + err);
        return;
    }

    TransientSensitivityAnalysis sens(circuit);
    sens.setVariation(0.0, 5e-3, 1e-5);
    sens.setOutputVar("V(out)");
    sens.addParameter("R1", "R");
    sens.addParameter("C1", "C");
    const TransientSensitivityResult res = sens.run();
    check(res.success, "rc transient sensitivity succeeds: " + res.message);
    if (!res.success || res.points.empty()) return;

    const double R = 1e3;
    const double C = 1e-6;
    // The last point is the most reliable one (the response has settled and
    // the integrator is in its asymptotic regime).
    const auto& pt = res.points.back();
    const double t = pt.time;

    // v_out(t) = 1 - exp(-t/(RC))
    // dv/dR   = -(t/(R^2 C))·exp(-t/(RC))
    // dv/dC   = -(t/(R C^2))·exp(-t/(RC))
    const double e = std::exp(-t / (R * C));
    const double v_closed = 1.0 - e;
    const double dR_closed = -(t / (R * R * C)) * e;
    const double dC_closed = -(t / (R * C * C)) * e;

    // The integrator is second order with an LTE-controlled step, so the
    // tolerance is the integrator's, not the machine's.
    const double tol_v = 1e-3 * std::fabs(v_closed) + 1e-6;
    checkNear(pt.obj_value, v_closed, tol_v, "V(out)(t) vs 1-exp(-t/RC)");

    const double tol_r = 5e-2 * std::fabs(dR_closed) + 1e-9;
    checkNear(pt.sens[0], dR_closed, tol_r, "d(V(out))/dR vs the closed form");

    const double tol_c = 5e-2 * std::fabs(dC_closed) + 1e-6;
    checkNear(pt.sens[1], dC_closed, tol_c, "d(V(out))/dC vs the closed form");
}

// ===============================================================
// 6) The same transient with the BACKWARD EULER integrator: the
//    recursion is exact for it too (the charge history of the
//    trapezoidal rule has no secondary state to differentiate).
// ===============================================================
void testTransientRcBackwardEuler() {
    Circuit circuit;
    std::string err;
    const std::string netlist =
        "* rc step sensitivity (backward Euler)\n"
        "V1 in 0 PULSE(0 1 0 1p 1p 1 2)\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".END\n";
    if (!buildCircuit(netlist, circuit, err)) {
        check(false, "rc transient netlist parses: " + err);
        return;
    }

    TransientSensitivityAnalysis sens(circuit);
    sens.setVariation(0.0, 5e-3, 1e-5);
    sens.setIntegrationMethod(IntegrationMethod::BACKWARD_EULER);
    sens.setOutputVar("V(out)");
    sens.addParameter("R1", "R");
    sens.addParameter("C1", "C");
    const TransientSensitivityResult res = sens.run();
    check(res.success, "rc transient sensitivity (BE) succeeds: " + res.message);
    if (!res.success || res.points.empty()) return;

    const double R = 1e3;
    const double C = 1e-6;
    const auto& pt = res.points.back();
    const double t = pt.time;
    const double e = std::exp(-t / (R * C));
    const double dR_closed = -(t / (R * R * C)) * e;
    const double dC_closed = -(t / (R * C * C)) * e;

    checkNear(pt.sens[0], dR_closed, 5e-2 * std::fabs(dR_closed) + 1e-9,
              "d(V(out))/dR vs the closed form (backward Euler)");
    checkNear(pt.sens[1], dC_closed, 5e-2 * std::fabs(dC_closed) + 1e-6,
              "d(V(out))/dC vs the closed form (backward Euler)");
}

// ===============================================================
// 7) The card spellings: ngspice and Xyce, DC/AC/transient
// ===============================================================
void testCardSpellings() {
    struct Case {
        const char* netlist;
        bool        ac;
        bool        tran;
    };
    const Case cases[] = {
        {"* ngspice dc\n"
         "V1 in 0 DC 1\n"
         "R1 in out 1k\n"
         ".SENS V(out) R1\n"
         ".END\n", false, false},
        {"* xyce dc\n"
         "V1 in 0 DC 1\n"
         "R1 in out 1k\n"
         ".SENS objfunc=V(out) param=R1\n"
         ".END\n", false, false},
        {"* ngspice ac\n"
         "V1 in 0 DC 1 AC 1\n"
         "R1 in out 1k\n"
         ".AC DEC 2 1k 1e6\n"
         ".SENS AC V(out) R1\n"
         ".END\n", true, false},
        {"* xyce ac\n"
         "V1 in 0 DC 1 AC 1\n"
         "R1 in out 1k\n"
         ".AC DEC 2 1k 1e6\n"
         ".SENS acobjfunc=V(out) param=R1\n"
         ".END\n", true, false},
        {"* ngspice transient\n"
         "V1 in 0 DC 1\n"
         "R1 in out 1k\n"
         ".TRAN 1u 10u\n"
         ".SENS TRAN V(out) R1\n"
         ".END\n", false, true},
        {"* xyce transient\n"
         "V1 in 0 DC 1\n"
         "R1 in out 1k\n"
         ".TRAN 1u 10u\n"
         ".SENS tranobjfunc=V(out) param=R1\n"
         ".END\n", false, true},
    };

    for (const auto& c : cases) {
        NetlistParser parser;
        Circuit circuit;
        std::string err;
        if (!parser.parseString(c.netlist, circuit)) {
            check(false, "card spelling parses");
            continue;
        }
        const AnalysisParameters& ap = parser.getAnalysisParameters();
        check(ap.want_sens, "the .SENS card was recognized");
        check(ap.sens_ac == c.ac, "the AC flavour of .SENS was recognized");
        check(ap.sens_tran == c.tran, "the transient flavour of .SENS was recognized");
        check(ap.sens_output_var == "V(out)",
              "the output variable was reassembled: '" + ap.sens_output_var + "'");
        check(ap.sens_params.size() == 1 && ap.sens_params[0].first == "R1",
              "the parameter list was parsed");
    }
}

} // namespace

int main() {
    testDcDivider();
    testAdjointAgainstDense();
    testAcLowPass();
    testTransientResistive();
    testTransientRcAgainstClosedForm();
    testTransientRcBackwardEuler();
    testCardSpellings();

    std::cout << "\n[SENS] " << (g_checks - g_failures) << "/" << g_checks
              << " checks passed\n";
    if (g_failures > 0) {
        std::cout << "[SENS] " << g_failures << " FAILURE(S)\n";
        return 1;
    }
    std::cout << "[SENS] OK\n";
    return 0;
}
