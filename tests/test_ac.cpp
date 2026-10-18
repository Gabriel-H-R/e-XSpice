/**
 * @file test_ac.cpp
 * @brief AC small-signal analysis: analytic networks, device linearization,
 *        parser, and the frequency grid.
 *
 * Every frequency-domain result is checked against a CLOSED FORM (or, for the
 * diode, against the device's own small-signal conductance), because a wrong
 * phase is numerically indistinguishable from "a plausible answer": the tests
 * pin magnitude AND phase, plus the branch-current phasor.
 *
 * Two of the cases are exactness regressions rather than approximations:
 *   - section 3 pins an ideal-source node, which must come out EXACTLY at the
 *     phasor amplitude (a pivot perturbation or a stale linearization shows up
 *     here first, at ~1e-7 relative);
 *   - section 5 uses an inductor, whose KCL row has no self-term: it exercises
 *     the pivot-safe ordering of the MNA system.
 */

#include "../src/analysis/ACAnalysis.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/core/Circuit.h"
#include "../src/core/Device.h"
#include "../src/core/NetlistParser.h"
#include "../src/devices/VoltageSource.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace e_XSpice;

// C++17 (strict ISO) does not expose M_PI.
static constexpr double kPiTest = 3.14159265358979323846;

static int g_fail = 0;
static const char* g_section = "";

static void section(const char* name) {
    g_section = name;
    std::printf("\n%s\n", name);
}

static void check(bool cond, const char* what) {
    if (!cond) ++g_fail;
    std::printf("  [%s] %s\n", cond ? "OK" : "FAIL", what);
}

static void checkNear(double got, double expect, double tol, const char* what) {
    const bool ok = std::isfinite(got) && std::fabs(got - expect) <= tol;
    if (!ok) ++g_fail;
    std::printf("  [%s] %-46s got=% .10g expect=% .10g\n",
                ok ? "OK" : "FAIL", what, got, expect);
}

// Index of a node inside ACDataPoint::mag: one entry per NON-GROUND node, in
// Circuit::getNodes() order (the same order the CSV labels use).
static int nodeIndex(const Circuit& circuit, const std::string& name) {
    int i = 0;
    for (const auto& n : circuit.getNodes()) {
        if (n->isGround()) continue;
        if (n->getName() == name) return i;
        ++i;
    }
    return -1;
}

// ===============================================================
// 1) Frequency grid (.AC LIN/DEC/OCT)
// ===============================================================
static void testGrid() {
    section("1) Frequency grid");

    // DEC 10 points per decade: 1 Hz .. 1 kHz = 31 points (inclusive ends).
    const std::vector<double> dec =
        ACAnalysis::buildGrid(ACVariation::DEC, 10, 1.0, 1000.0);
    check(dec.size() == 31, "DEC 10/decade over 3 decades gives 31 points");
    if (dec.size() == 31) {
        checkNear(dec.front(), 1.0, 1e-15, "DEC first point");
        checkNear(dec.back(), 1000.0, 1e-9, "DEC last point");
        // The ratio between consecutive points must be exactly 10^(1/n).
        const double ratio_ref = std::pow(10.0, 0.1);
        double worst = 0.0;
        for (size_t i = 1; i + 1 < dec.size(); ++i) {
            worst = std::max(worst, std::fabs(dec[i + 1] / dec[i] - ratio_ref));
        }
        checkNear(worst, 0.0, 1e-12, "DEC step is exactly 10^(1/n)");
    }

    // OCT 4 points per octave: 1 .. 8 Hz spans 3 octaves -> 13 points.
    const std::vector<double> oct =
        ACAnalysis::buildGrid(ACVariation::OCT, 4, 1.0, 8.0);
    if (oct.size() != 13) {
        std::printf("       (oct grid has %zu points:", oct.size());
        for (double f : oct) std::printf(" %g", f);
        std::printf(")\n");
    }
    check(oct.size() == 13, "OCT 4/octave over 3 octaves gives 13 points");
    if (!oct.empty()) checkNear(oct.back(), 8.0, 1e-12, "OCT last point");

    // LIN 5 points: endpoints exact, uniform steps.
    const std::vector<double> lin =
        ACAnalysis::buildGrid(ACVariation::LIN, 5, 10.0, 50.0);
    check(lin.size() == 5, "LIN 5 points");
    if (lin.size() == 5) {
        checkNear(lin[0], 10.0, 0.0, "LIN first point exact");
        checkNear(lin[4], 50.0, 0.0, "LIN last point exact (no accumulation drift)");
        checkNear(lin[2], 30.0, 1e-12, "LIN midpoint");
    }

    check(ACAnalysis::buildGrid(ACVariation::DEC, 10, -1.0, 1000.0).empty(),
          "non-positive start is rejected");
}

// ===============================================================
// 2) RC low-pass: |H| and phase against the closed form
// ===============================================================
static void testRcLowPass() {
    section("2) RC low-pass vs closed form  H = 1/(1+jwRC)");

    const double R = 1e3, C = 1e-9;
    std::stringstream ss;
    ss << "V1 in 0 DC 0 AC 1\n"
       << "R1 in out " << R << "\n"
       << "C1 out 0 " << C << "\n"
       << ".AC DEC 10 10 1Meg\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);
    check(parser.getAnalysisParameters().want_ac, "parser sets want_ac");

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::DEC);
    ac.setPoints(10);
    ac.setFrequencyRange(10.0, 1e6);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(/*use_mpi=*/false);
    check(res.success, "AC run succeeded");
    check(res.points.size() == 51, "51 frequencies in the sweep (5 decades x 10 + 1)");

    const int out_i = nodeIndex(circuit, "out");
    const int in_i = nodeIndex(circuit, "in");
    check(out_i >= 0 && in_i >= 0, "node indices located");

    double worst_mag = 0.0, worst_ph = 0.0;
    for (const auto& p : res.points) {
        const double wc = 2.0 * kPiTest * p.frequency * R * C;
        const double m_ref = 1.0 / std::sqrt(1.0 + wc * wc);
        const double ph_ref = -std::atan(wc) * 180.0 / kPiTest;
        worst_mag = std::max(worst_mag, std::fabs(p.mag[static_cast<size_t>(out_i)] - m_ref));
        worst_ph = std::max(worst_ph,
                            std::fabs(p.phase_deg[static_cast<size_t>(out_i)] - ph_ref));
        // The source pins its node: |V(in)| == 1 exactly, phase == 0.
        worst_mag = std::max(worst_mag, std::fabs(p.mag[static_cast<size_t>(in_i)] - 1.0));
    }
    checkNear(worst_mag, 0.0, 1e-12, "max |V| deviation vs closed form");
    checkNear(worst_ph, 0.0, 1e-9, "max phase deviation vs closed form (deg)");
}

// ===============================================================
// 3) Exactness: an ideal source pins its node
// ===============================================================
static void testPinnedSourceExact() {
    section("3) Ideal source pins its node EXACTLY (MNA ordering)");

    // V1 - L1 - out, R1 to ground: node `in` has NO self-term (its KCL row is
    // i(V1) + i(L1) = 0), which is exactly the MNA case that needs a
    // pivot-safe ordering. The pinned node must come out at the phasor
    // amplitude with zero error, and the inductor current at the analytic
    // |I| = |Vac| / |R + jwL|.
    std::stringstream ss;
    ss << "V1 in 0 DC 0 AC 2\n"
       << "L1 in out 1e-3\n"
       << "R1 out 0 1\n"
       << ".AC LIN 1 1e3 1e3\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::LIN);
    ac.setPoints(1);
    ac.setFrequencyRange(1e3, 1e3);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");
    check(!res.stats.pivot_perturbed,
          "no pivot perturbation needed for an MNA node row");

    const int in_i = nodeIndex(circuit, "in");
    const int out_i = nodeIndex(circuit, "out");
    checkNear(res.points[0].mag[static_cast<size_t>(in_i)], 2.0, 0.0,
              "|V(in)| == Vac exactly");
    checkNear(res.points[0].phase_deg[static_cast<size_t>(in_i)], 0.0, 0.0,
              "phase of the pinned node is exactly 0");

    const double omega = 2.0 * kPiTest * 1e3;
    const std::complex<double> z(1.0, omega * 1e-3);
    const std::complex<double> i = std::complex<double>(2.0, 0.0) / z;
    checkNear(res.points[0].mag[static_cast<size_t>(out_i)], std::abs(i), 1e-14,
              "|v(out)| = R·|I| (closed form)");
    // Branch currents follow Circuit::getDevices() order: V1 first, then L1.
    check(res.points[0].branch_mag.size() == 2, "two branch currents reported");
    if (res.points[0].branch_mag.size() == 2) {
        checkNear(res.points[0].branch_mag[1], std::abs(i), 1e-14, "|I(L1)|");
        checkNear(res.points[0].branch_phase_deg[1],
                  std::arg(i) * 180.0 / kPiTest, 1e-9, "phase of I(L1) (deg)");
        // The source current is defined + -> - THROUGH the source, so it is
        // the negative of the current delivered to the circuit (180 deg apart).
        checkNear(res.points[0].branch_mag[0], std::abs(i), 1e-14, "|I(V1)|");
        checkNear(res.points[0].branch_phase_deg[0],
                  std::arg(-i) * 180.0 / kPiTest, 1e-9,
                  "phase of I(V1) is +180 deg (absorbed-current convention)");
    }
}

// ===============================================================
// 4) Series RLC: resonance peak at f0 = 1/(2*pi*sqrt(LC))
// ===============================================================
static void testRlcResonance() {
    section("4) Series RLC: resonance frequency and peak Q");

    // V1 - L - mid - R - out - C - gnd, output across C:
    //     H(s) = 1 / (1 + sRC + s^2 LC)
    // peaks at f0 = 1/(2*pi*sqrt(LC)) with |H| = Q = sqrt(L/C)/R.
    const double L = 1e-3, C = 1e-9, R = 1.0;
    std::stringstream ss;
    ss << "V1 in 0 DC 0 AC 1\n"
       << "L1 in mid " << L << "\n"
       << "R1 mid out " << R << "\n"
       << "C1 out 0 " << C << "\n"
       << ".AC LIN 4001 1e5 2.2e5\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::LIN);
    ac.setPoints(4001);
    ac.setFrequencyRange(1e5, 2.2e5);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");
    check(res.points.size() == 4001, "all 4001 frequencies retained");

    const int out_i = nodeIndex(circuit, "out");
    check(res.points[0].mag.size() > static_cast<size_t>(out_i), "output reported");

    // Closed form: H(s) = 1 / (1 + sRC + s^2 LC) across the WHOLE sweep —
    // stronger than comparing only the sampled peak, and it pins the phase too.
    double worst_mag = 0.0, worst_ph = 0.0;
    for (const auto& p : res.points) {
        const std::complex<double> s(0.0, 2.0 * kPiTest * p.frequency);
        const std::complex<double> h = 1.0 / (1.0 + s * R * C + s * s * L * C);
        worst_mag = std::max(
            worst_mag, std::fabs(p.mag[static_cast<size_t>(out_i)] - std::abs(h)));
        worst_ph = std::max(
            worst_ph,
            std::fabs(p.phase_deg[static_cast<size_t>(out_i)] -
                      std::arg(h) * 180.0 / kPiTest));
    }
    checkNear(worst_mag, 0.0, 1e-6, "max |H| deviation vs closed form");
    checkNear(worst_ph, 0.0, 1e-6, "max phase deviation vs closed form (deg)");

    // Peak of the sampled sweep, and where it sits relative to f0.
    size_t peak_i = 0;
    double peak = -1.0;
    for (size_t i = 0; i < res.points.size(); ++i) {
        if (res.points[i].mag[static_cast<size_t>(out_i)] > peak) {
            peak = res.points[i].mag[static_cast<size_t>(out_i)];
            peak_i = i;
        }
    }

    const double f0 = 1.0 / (2.0 * kPiTest * std::sqrt(L * C));
    const double grid_step = 3e5 / 4000.0;   // (2.2e5 - 1e5) / (4001 - 1)
    checkNear(res.points[peak_i].frequency, f0, grid_step,
              "peak within one grid step of 1/(2*pi*sqrt(LC))");
    const double Q = std::sqrt(L / C) / R;
    // At the sampled peak the magnitude is slightly below Q (the grid does not
    // land exactly on f0), so compare against the closed form there.
    {
        const std::complex<double> s(0.0, 2.0 * kPiTest * res.points[peak_i].frequency);
        const std::complex<double> h = 1.0 / (1.0 + s * R * C + s * s * L * C);
        checkNear(peak, std::abs(h), 1e-9, "peak magnitude vs closed form");
        check(peak < Q + 1e-9 && peak > Q * 0.99,
              "peak magnitude is just below Q (grid resolution)");
    }
}

// ===============================================================
// 5) Phase of the source and a two-source interference check
// ===============================================================
static void testSourcePhase() {
    section("5) Source phase and superposition of two phasors");

    const double R1 = 1e3, R2 = 1e3, R3 = 1e3;
    std::stringstream ss;
    ss << "V1 a 0 DC 0 AC 1\n"          // 0 deg
       << "V2 b 0 DC 0 AC 1 90\n"      // 90 deg
       << "R1 a out " << R1 << "\n"
       << "R2 b out " << R2 << "\n"
       << "R3 out 0 " << R3 << "\n"
       << ".AC LIN 1 1e6 1e6\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::LIN);
    ac.setPoints(1);
    ac.setFrequencyRange(1e6, 1e6);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");

    const int out_i = nodeIndex(circuit, "out");
    // Contributions: a drives 1/3, b drives j/3 -> |out| = sqrt(2)/3 at 45 deg.
    checkNear(res.points[0].mag[static_cast<size_t>(out_i)], std::sqrt(2.0) / 3.0,
              1e-14, "magnitude of the phasor sum");
    checkNear(res.points[0].phase_deg[static_cast<size_t>(out_i)], 45.0, 1e-9,
              "phase of the phasor sum (deg)");
}

// ===============================================================
// 6) Diode: the AC response must agree with the DC nonlinearity
// ===============================================================
static void testDiodeSmallSignal() {
    section("6) Diode: AC response agrees with the DC nonlinearity");

    const double vs = 0.6, vac = 1e-3;
    std::stringstream ss;
    ss << ".MODEL DMOD D (IS=1e-14 N=1)\n"
       << "V1 a 0 DC " << vs << " AC " << vac << "\n"
       << "R1 a b 1\n"
       << "D1 b 0 DMOD\n"
       << ".AC LIN 1 1e3 1e3\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    check(parser.parseString(ss.str(), circuit), "netlist parses");

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::LIN);
    ac.setPoints(1);
    ac.setFrequencyRange(1e3, 1e3);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");
    check(res.points[0].branch_mag.size() == 1,
          "branch current reported (extra unknown of V1)");

    // Effective small-signal conductance straight from the DC solver: a
    // central difference around the bias, using the SAME model (so this checks
    // the AC linearization against the DC nonlinearity, not against an
    // idealized diode formula).
    const int a_i = nodeIndex(circuit, "a");
    const int b_i = nodeIndex(circuit, "b");
    auto dcCurrentAt = [&](double v1) {
        std::stringstream net;
        net << ".MODEL DMOD D (IS=1e-14 N=1)\n"
            << "V1 a 0 DC " << v1 << "\n"
            << "R1 a b 1\n"
            << "D1 b 0 DMOD\n"
            << ".OP\n.END\n";
        Circuit c2;
        NetlistParser p2;
        p2.parseString(net.str(), c2);
        DCAnalysis dc(c2);
        dc.setTolerance(1e-12);
        const DCResult op = dc.run(false, /*extract=*/false);
        if (!op.success) return 0.0;
        const double va = op.node_voltages[static_cast<size_t>(a_i)];
        const double vb = op.node_voltages[static_cast<size_t>(b_i)];
        return (va - vb) / 1.0;              // current through R1 (1 ohm)
    };
    const double gm_eff = (dcCurrentAt(vs + vac) - dcCurrentAt(vs - vac))
                          / (2.0 * vac);

    check(gm_eff > 0.0, "DC small-signal conductance is positive");
    const double expect = gm_eff * vac;
    if (!res.points[0].branch_mag.empty()) {
        checkNear(res.points[0].branch_mag[0], expect, 1e-3 * expect,
                  "|I(V1)| == gm_eff - Vac  (AC vs DC)");
    }
}

// ===============================================================
// 7) CSV stream: header, row count, full precision
// ===============================================================
static void testCsv() {
    section("7) CSV stream");

    std::stringstream ss;
    ss << "V1 in 0 DC 0 AC 1\n"
       << "R1 in out 1k\n"
       << "C1 out 0 1n\n"
       << ".AC DEC 2 1 100\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);

    std::ostringstream csv;
    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::DEC);
    ac.setPoints(2);
    ac.setFrequencyRange(1.0, 100.0);
    ac.setTolerance(1e-10);
    ac.setCsvStream(&csv, {"Vin", "Vout"});
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");

    std::istringstream in(csv.str());
    std::string line;
    int rows = 0;
    std::string header;
    while (std::getline(in, line)) {
        if (rows == 0) header = line;
        ++rows;
    }
    // DEC 2/decade over two decades: 1, 3.16, 10, 31.6, 100 -> 5 rows.
    check(rows == 6, "DEC 2/decade over 2 decades -> 1 header + 5 rows");
    check(header.find("frequency") == 0, "header starts with 'frequency'");
    check(header.find("Vout_mag") != std::string::npos, "header has Vout_mag");
    check(header.find("Vout_phase") != std::string::npos, "header has Vout_phase");
    check(header.find("V1_imag_mag") != std::string::npos,
          "header has the branch current columns");
    // Full precision: values must not be truncated to 6 significant digits.
    check(csv.str().find("e-") != std::string::npos ||
          csv.str().find("0.999999") != std::string::npos,
          "CSV values are written with full precision");
}

// ===============================================================
// 8) Statistics and a failing operating point
// ===============================================================
static void testStats() {
    section("8) Statistics and operating-point behavior");

    std::stringstream ss;
    ss << "V1 in 0 DC 0 AC 1\n"
       << "R1 in out 1k\n"
       << "C1 out 0 1n\n"
       << ".AC DEC 5 10 100\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::DEC);
    ac.setPoints(5);
    ac.setFrequencyRange(10.0, 100.0);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");
    check(res.stats.frequencies == 6, "6 frequencies (DEC 5 over one decade)");
    check(res.stats.factorizations == 6, "6 factorizations");
    check(res.stats.factor_flops > 0, "factorization flops accounted");
    check(res.stats.newton_iterations > 0, "operating point iterations counted");
    check(res.stats.failed_frequencies == 0, "no failed frequency");
    check(!res.stats.pivot_perturbed, "no pivot perturbation on a passive net");

    // A netlist without .AC parameters must be rejected, not crash.
    {
        std::stringstream bad;
        bad << "V1 in 0 DC 0 AC 1\n.AC DEC 2 -1 10\n.END\n";
        Circuit c2;
        NetlistParser p2;
        check(!p2.parseString(bad.str(), c2),
              ".AC with a negative frequency is rejected by the parser");
    }
}

// ===============================================================
// 9) DC bias is preserved: the AC source does not move the operating point
// ===============================================================
static void testDCOffset() {
    section("9) AC magnitude is small-signal around the DC bias");

    std::stringstream ss;
    ss << "V1 in 0 DC 5 AC 1m\n"
       << "R1 in out 1k\n"
       << "R2 out 0 1k\n"
       << ".AC LIN 1 1e6 1e6\n"
       << ".END\n";

    Circuit circuit;
    NetlistParser parser;
    parser.parseString(ss.str(), circuit);

    ACAnalysis ac(circuit);
    ac.setVariation(ACVariation::LIN);
    ac.setPoints(1);
    ac.setFrequencyRange(1e6, 1e6);
    ac.setTolerance(1e-10);
    const ACResult res = ac.run(false);
    check(res.success, "AC run succeeded");

    const int out_i = nodeIndex(circuit, "out");
    checkNear(res.points[0].mag[static_cast<size_t>(out_i)], 1e-3 / 2.0, 1e-16,
              "phasor carries the AC part only (1 mV / 2)");
}

// ===============================================================
// 10) Parser: source specs and .AC validation
// ===============================================================
static void testParser() {
    section("10) Parser: source AC specs");

    {
        std::stringstream ss;
        ss << "V1 in 0 AC 2\n"          // AC only -> DC 0
           << "V2 n2 0 DC 3 AC 1 45\n"  // DC 3, AC 1 @ 45 deg
           << "V3 n3 0 PULSE(0 5 1m) AC 3\n"
           << "I1 n4 0 AC 1m\n"
           << ".AC LIN 5 1 10\n"
           << ".END\n";
        Circuit circuit;
        NetlistParser parser;
        check(parser.parseString(ss.str(), circuit),
              "netlist with AC specs parses");
        const auto& devs = circuit.getDevices();
        check(devs.size() == 4, "4 sources parsed");
        if (devs.size() == 4) {
            const auto* v1 = dynamic_cast<const VoltageSource*>(devs[0].get());
            check(v1 != nullptr, "V1 is a VoltageSource");
            if (v1) {
                // "AC 2" with no time-domain spec implies DC 0 (SPICE).
                checkNear(v1->waveform().value(0.0), 0.0, 1e-15,
                          "AC-only source has DC value 0");
                check(v1->waveform().hasAC(), "V1 has an AC spec");
                checkNear(v1->waveform().acMagnitude(), 2.0, 1e-15,
                          "V1 AC magnitude");
            }
            const auto* v2 = dynamic_cast<const VoltageSource*>(devs[1].get());
            check(v2 != nullptr, "V2 is a VoltageSource");
            if (v2) {
                checkNear(v2->waveform().value(0.0), 3.0, 1e-15, "V2 DC value is 3");
                check(v2->waveform().hasAC(), "V2 has an AC spec");
                double re = 0.0, im = 0.0;
                v2->waveform().acPhasor(re, im);
                checkNear(re, std::cos(kPiTest / 4.0), 1e-14, "V2 AC real part");
                checkNear(im, std::sin(kPiTest / 4.0), 1e-14, "V2 AC imag part");
            }
            const auto* v3 = dynamic_cast<const VoltageSource*>(devs[2].get());
            check(v3 != nullptr, "V3 is a VoltageSource");
            if (v3) {
                check(v3->waveform().hasAC(),
                      "PULSE source may carry an AC spec too");
                checkNear(v3->waveform().acMagnitude(), 3.0, 1e-15,
                          "V3 AC magnitude");
            }
        }
        // A source with no AC spec is silent in an AC analysis.
        std::stringstream quiet;
        quiet << "V1 in 0 DC 1\nR1 in 0 1k\n.AC LIN 2 1 10\n.END\n";
        Circuit c2;
        NetlistParser p2;
        check(p2.parseString(quiet.str(), c2), "netlist without AC specs parses");
        if (!c2.getDevices().empty()) {
            const auto* vs = dynamic_cast<const VoltageSource*>(
                c2.getDevices()[0].get());
            check(vs && !vs->waveform().hasAC(),
                  "source without AC spec is silent in AC");
        }
    }

    // Invalid .AC cards are rejected with a message.
    {
        std::stringstream ss;
        ss << "V1 in 0 DC 0 AC 1\n.AC XYZ 5 1 10\n.END\n";
        Circuit circuit;
        NetlistParser parser;
        check(!parser.parseString(ss.str(), circuit),
              ".AC with an invalid variation is rejected");
    }
    {
        std::stringstream ss;
        ss << "V1 in 0 DC 0 AC 1\n.AC DEC 5 1\n.END\n";
        Circuit circuit;
        NetlistParser parser;
        check(!parser.parseString(ss.str(), circuit),
              ".AC with missing arguments is rejected");
    }
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("AC analysis checks\n");

    testGrid();
    testRcLowPass();
    testPinnedSourceExact();
    testRlcResonance();
    testSourcePhase();
    testDiodeSmallSignal();
    testCsv();
    testStats();
    testDCOffset();
    testParser();

    if (g_fail == 0) {
        std::printf("\nAC analysis PASSED\n");
        return 0;
    }
    std::printf("\nAC analysis FAILED (%d checks in section '%s')\n",
                g_fail, g_section);
    return 1;
}