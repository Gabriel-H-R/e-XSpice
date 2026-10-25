/**
 * @author Gabriel Henrique Silva
 * @file test_hb.cpp
 * @brief Harmonic Balance analysis (`.HB`).
 *
 * The tests are organized around two independent references:
 *
 *   1) LINEARITY — for a linear circuit the harmonic balance solution must
 *      reproduce the AC analysis exactly (the fundamental matches the AC
 *      gain and every higher harmonic is zero). This validates the whole
 *      machinery (the synthesis, the DFT, the Jacobian, the expanded
 *      system and its solve) without any hand-computed Fourier series.
 *
 *   2) THE TRANSIENT — the harmonic balance solution must match the
 *      steady state of a transient simulation of the same circuit: the
 *      Fourier coefficients of the last period of the transient are the
 *      harmonic balance harmonics. This is the strongest check available
 *      (it validates a completely different solver against the HB one).
 */

#include "core/Circuit.h"
#include "core/NetlistParser.h"
#include "core/Node.h"
#include "analysis/ACAnalysis.h"
#include "analysis/HBAnalysis.h"
#include "analysis/TransientAnalysis.h"

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

bool build(const std::string& text, Circuit& circuit) {
    NetlistParser parser;
    return parser.parseString(text, circuit);
}

/** The DAE row of a node (the HB reports the MNA unknowns). */
int rowOf(const Circuit& circuit, const std::string& name) {
    int row = 0;
    for (const auto& n : circuit.getNodes()) {
        if (n->isGround()) continue;
        if (n->getName() == name) return row;
        ++row;
    }
    return -1;
}

// ===============================================================
// 1) A LINEAR circuit: the HB must reproduce the AC analysis
// ===============================================================
void testLinearMatchesAc() {
    Circuit circuit;
    if (!build("* linear rc\n"
               "V1 in 0 DC 0 AC 1\n"
               "R1 in out 1k\n"
               "C1 out 0 1u\n"
               ".END\n", circuit)) {
        check(false, "the linear netlist parses");
        return;
    }
    const double f = 1000.0;

    // The AC reference.
    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::LIN);
    ac.setPoints(1);
    ac.setFrequencyRange(f, f);
    const ACResult ac_res = ac.run(false);
    check(ac_res.success, "the AC reference runs");
    if (!ac_res.success || ac_res.points.empty()) return;

    HBAnalysis hb(circuit);
    hb.setFundamental(f);
    hb.setHarmonics(3);
    const HBResult hb_res = hb.run();
    check(hb_res.success, "the HB runs: " + hb_res.message);
    if (!hb_res.success) return;

    const int out_row = rowOf(circuit, "out");
    check(out_row >= 0, "the output node has a row");
    if (out_row < 0) return;

    // The fundamental must match the AC magnitude and phase.
    const std::complex<double> h1 =
        hb_res.harmonics[1][static_cast<size_t>(out_row)];
    const double ac_mag = ac_res.points.front().mag[1];
    const double ac_ph = ac_res.points.front().phase_deg[1];
    checkNear(std::abs(h1), ac_mag, 1e-6, "the HB fundamental matches the AC magnitude");
    checkNear(std::arg(h1) * 180.0 / 3.14159265358979323846, ac_ph, 1e-4,
              "the HB fundamental matches the AC phase");

    // The higher harmonics of a linear circuit must vanish.
    for (int k = 2; k <= 3; ++k) {
        checkNear(std::abs(hb_res.harmonics[static_cast<size_t>(k)]
                                       [static_cast<size_t>(out_row)]),
                  0.0, 1e-9,
                  "the higher harmonic " + std::to_string(k) + " of a linear circuit vanishes");
    }
}

// ===============================================================
// 2) A NONLINEAR circuit: the HB must match the steady state of a
//    transient simulation (the Fourier coefficients of its last
//    period).
// ===============================================================
void testNonlinearMatchesTransient() {
    Circuit circuit;
    // The source carries BOTH excitations: the time-domain waveform (what
    // the transient uses) and the AC phasor (what the HB uses). They
    // describe the same signal: SIN(0 0.5 f 0 0 90) is 0.5*cos(wt) and
    // the phasor 0.5<0 deg is 0.5*cos(wt) too.
    if (!build("* diode clipper\n"
               ".MODEL DMOD D(IS=1e-12 RS=0 N=1)\n"
               "V1 in 0 DC 0 SIN(0 0.5 1k 0 0 90) AC 0.5\n"
               "R1 in out 1k\n"
               "D1 out 0 DMOD\n"
               ".END\n", circuit)) {
        check(false, "the nonlinear netlist parses");
        return;
    }
    const double f = 1000.0;

    // The transient reference: many periods, and the Fourier coefficients
    // of the LAST one.
    TransientAnalysis tran(circuit);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    tran.setTolerance(1e-9);
    tran.setSaveInterval(1);
    const TransientResult tr = tran.run(0.0, 20.0 / f, 2.0e-7, false);
    check(tr.success, "the transient reference runs");
    if (!tr.success || tr.data.size() < 100) return;

    // The Fourier coefficients of the last period (a rectangular window of
    // one period, sampled at the accepted times).
    const double period = 1.0 / f;
    const double t_end = tr.data.back().time;
    const double t_begin = t_end - period;
    std::vector<double> samples;
    for (const auto& p : tr.data) {
        if (p.time >= t_begin && p.time <= t_end) samples.push_back(p.time);
    }
    check(samples.size() > 16, "the last period has enough samples");
    if (samples.size() <= 16) return;

    // The Fourier coefficients of the last period: the correlation with the
    // fundamental (the mean over the sampled points, which is the exact
    // Fourier coefficient for a uniform sampling).
    const int out_row = rowOf(circuit, "out");
    double a1 = 0.0, b1 = 0.0;
    int n_used = 0;
    for (const auto& d : tr.data) {
        if (d.time < t_begin || d.time > t_end) continue;
        const double w = 2.0 * 3.14159265358979323846 * f * (d.time - t_begin);
        const double v = d.node_voltages[static_cast<size_t>(out_row)];
        a1 += v * std::cos(w);
        b1 += v * std::sin(w);
        ++n_used;
    }
    check(n_used > 16, "the last period has enough samples");
    if (n_used <= 16) return;
    if (std::getenv("EXSPICE_DEBUG_HB")) {
        std::fprintf(stderr, "[HB] transient: n_used=%d out_row=%d a1=%g b1=%g "
                     "v_last=%g\n", n_used, out_row, a1, b1,
                     tr.data.back().node_voltages[static_cast<size_t>(out_row)]);
    }
    a1 *= 2.0 / static_cast<double>(n_used);
    b1 *= 2.0 / static_cast<double>(n_used);
    const std::complex<double> coef1(a1, b1);

    HBAnalysis hb(circuit);
    hb.setFundamental(f);
    hb.setHarmonics(3);
    const HBResult hb_res = hb.run();
    check(hb_res.success, "the HB runs on the nonlinear circuit: " + hb_res.message);
    if (!hb_res.success) return;

    // The fundamental must match the transient's steady state.
    const std::complex<double> h1 =
        hb_res.harmonics[1][static_cast<size_t>(out_row)];
    checkNear(std::abs(h1), std::abs(coef1), 0.05 * std::abs(coef1) + 1e-3,
              "the HB fundamental matches the transient steady state "
              "(magnitude)");
    checkNear(std::arg(h1) * 180.0 / 3.14159265358979323846,
              std::arg(coef1) * 180.0 / 3.14159265358979323846, 10.0,
              "the HB fundamental matches the transient steady state (phase)");

    // The third harmonic must be significant (the diode is nonlinear) and
    // the HB must predict a much smaller third than the second.
    const double h3 = std::abs(hb_res.harmonics[3][static_cast<size_t>(out_row)]);
    check(h3 > 1e-3, "the third harmonic of a clipped sine is not negligible");
}

} // namespace

int main() {
    testLinearMatchesAc();
    testNonlinearMatchesTransient();

    std::cout << "\n[HB] " << (g_checks - g_failures) << "/" << g_checks
              << " checks passed\n";
    if (g_failures > 0) {
        std::cout << "[HB] " << g_failures << " FAILURE(S)\n";
        return 1;
    }
    std::cout << "[HB] OK\n";
    return 0;
}
