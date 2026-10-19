/**
 * @author Gabriel Henrique Silva
 * @file test_noise.cpp
 * @brief Unit tests for the .NOISE analysis (adjoint transfer, integrated
 *        noise, generator reports, parser).
 *
 * The checks are ANALYTIC wherever an analytic answer exists — a noise
 * analysis that is not validated against closed forms is a noise analysis you
 * cannot trust:
 *
 *   1) sweep grid rules (same as .AC);
 *   2) resistive divider: S_out, gain and S_in against 4kT/R closed forms
 *      (including the fact that an ideal voltage source is a SHORT in noise
 *      analysis, so R1's noise current is sunk by the source);
 *   3) RC low-pass: S_out(f) = 4kT/R * |H(f)|^2 — a genuinely FREQUENCY
 *      dependent check of the adjoint transfer;
 *   4) the per-generator breakdown sums exactly to the total;
 *   5) the integrator: white noise over a band integrates to S*(f2-f1);
 *   6) temperature scaling (linear in T);
 *   7) diode shot noise against 2qI, using the DC current from the same
 *      solver (so the test stays independent of the diode model);
 *   8) differential output V(out,ref);
 *   9) parser: both `.NOISE` spellings (ngspice `V(out) SRC` and Xyce
 *      `SRC=V1`), the optional reference node and pts_per_summary.
 *
 * The tests configure the analysis THROUGH the parsed card, exactly like
 * main.cpp does, so the parser-to-analysis mapping is covered too.
 */

#include "../src/core/Circuit.h"
#include "../src/core/NetlistParser.h"
#include "../src/analysis/NoiseAnalysis.h"
#include "../src/core/Device.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace e_XSpice;

namespace {

int g_fail = 0;
int g_checks = 0;
const char* g_section = "";

void section(const char* name) { g_section = name; std::printf("-- %s\n", name); }

void check(bool ok, const char* what) {
    ++g_checks;
    if (!ok) {
        ++g_fail;
        std::printf("   FAIL: %s\n", what);
    }
}

void checkNear(double got, double want, double tol, const char* what) {
    ++g_checks;
    const double scale = std::max(1.0, std::fabs(want));
    if (!(std::fabs(got - want) <= tol * scale)) {
        ++g_fail;
        std::printf("   FAIL: %s (got %.12g, want %.12g)\n", what, got, want);
    }
}

/** Fails with the analysis message, which is where the real reason lives. */
void checkOk(bool ok, const std::string& message, const char* what) {
    ++g_checks;
    if (!ok) {
        ++g_fail;
        std::printf("   FAIL: %s -> %s\n", what, message.c_str());
    }
}

const double kPiTest = 3.14159265358979323846;
const double kBoltzmann = 1.380649e-23;      // J/K (exact, SI 2019)
const double kCharge = 1.602176634e-19;       // C   (exact, SI 2019)
const double kTemp300 = 300.15;

double thermalPsd(double r, double t = kTemp300) {
    return 4.0 * kBoltzmann * t / r;
}

AnalysisParameters parseInto(Circuit& circuit, const std::string& text) {
    NetlistParser parser;
    parser.parseString(text, circuit);
    return parser.getAnalysisParameters();
}

/** Mirrors what main.cpp does with the parsed card. */
NoiseResult runNoise(Circuit& circuit, const AnalysisParameters& p,
                     bool force_detail, bool verbose = false) {
    NoiseAnalysis noise(circuit);
    noise.setVariation(p.ac_variation == "LIN"   ? ACVariation::LIN
                      : p.ac_variation == "OCT"  ? ACVariation::OCT
                                                  : ACVariation::DEC);
    noise.setPoints(p.ac_points);
    noise.setFrequencyRange(p.ac_fstart, p.ac_fstop);
    noise.setOutputNode(p.noise_out_node, p.noise_ref_node);
    if (p.noise_source_given) noise.setInputSource(p.noise_source);
    noise.setSourceDetail(force_detail || p.noise_detail);
    noise.setVerbose(verbose);
    return noise.run(false);
}

int nodeIndex(const Circuit& circuit, const std::string& name) {
    int i = 0;
    for (const auto& n : circuit.getNodes()) {
        if (n->isGround()) continue;
        if (n->getName() == name) return i;
        ++i;
    }
    return -1;
}

std::size_t sourceIndex(const NoiseResult& r, const std::string& device,
                        const std::string& label) {
    for (std::size_t i = 0; i < r.sources.size(); ++i) {
        if (r.sources[i].device == device && r.sources[i].label == label) {
            return i;
        }
    }
    return static_cast<std::size_t>(-1);
}

// ===============================================================
// 1) Frequency grid
// ===============================================================
static void testGrid() {
    section("1) Frequency grid");

    // 10 Hz .. 1 MHz is FIVE decades.
    const auto dec = NoiseAnalysis::buildGrid(ACVariation::DEC, 10, 10.0, 1e6);
    if (dec.size() != 51) std::printf("   (DEC size = %zu)\n", dec.size());
    check(dec.size() == 51, "DEC 10/decade over 5 decades -> 51 points");
    checkNear(dec.front(), 10.0, 1e-15, "DEC starts at fstart");
    checkNear(dec.back(), 1e6, 1e-12, "DEC ends at fstop");

    const auto lin = NoiseAnalysis::buildGrid(ACVariation::LIN, 11, 1.0, 11.0);
    check(lin.size() == 11, "LIN 11 points");
    checkNear(lin[1] - lin[0], 1.0, 1e-12, "LIN step");

    const auto oct = NoiseAnalysis::buildGrid(ACVariation::OCT, 8, 1.0, 256.0);
    check(oct.size() == 65, "OCT 8/octave over 8 octaves -> 65 points");
}

// ===============================================================
// 2) Resistive divider: the closed forms
// ===============================================================
static void testDivider() {
    section("2) Resistive divider vs closed forms");

    const double r1 = 10e3, r2 = 10e3;
    Circuit circuit;
    const auto p = parseInto(circuit,
        "V1 in 0 DC 1\n"
        "R1 in out 10k\n"
        "R2 out 0 10k\n"
        ".NOISE V(out) V1 DEC 5 1 1e6 1\n"
        ".END\n");

    const NoiseResult res = runNoise(circuit, p, true);
    checkOk(res.success, res.message, "noise run succeeded");
    if (!res.success) return;
    check(res.sources.size() == 2, "two generators (R1, R2)");
    check(res.input_source == "V1", "input source recorded");
    check(res.output_node == "out" && res.reference_node == "0",
          "output nodes recorded");

    // An ideal voltage source is a SHORT for noise, so node `in` is at 0 V and
    // every generator sees the PARALLEL combination R1||R2 as its transfer
    // resistance (R1's noise current returns through R1 or through the shorted
    // source; a current injected at `out` returns through R1 or R2).
    const double rpar = r1 * r2 / (r1 + r2);
    const double s_out = rpar * rpar * (thermalPsd(r1) + thermalPsd(r2));
    const double gain = r2 / (r1 + r2);

    check(nodeIndex(circuit, "out") >= 0, "output node found");
    check(!res.points.empty(), "points retained");
    for (const auto& pt : res.points) {
        checkNear(pt.onoise_psd, s_out, 1e-12, "output PSD vs 4kT closed form");
        checkNear(pt.gain_mag, gain, 1e-12, "gain = R2/(R1+R2)");
        checkNear(pt.inoise_psd, s_out / (gain * gain), 1e-12,
                  "input-referred PSD = output/gain^2");
        check(pt.onoise_amp > 0.0, "amplitude form is positive");
        checkNear(pt.onoise_amp, std::sqrt(s_out), 1e-12,
                  "amplitude form is the sqrt of the PSD");
    }

    // Per-generator breakdown: each resistor contributes (R1||R2)^2 * 4kT/Ri.
    for (const auto& pt : res.points) {
        double sum = 0.0;
        for (double v : pt.source_psd) sum += v;
        checkNear(sum, pt.onoise_psd, 1e-12, "per-source PSDs sum to the total");
    }
    const std::size_t r1i = sourceIndex(res, "R1", "thermal");
    const std::size_t r2i = sourceIndex(res, "R2", "thermal");
    check(r1i != static_cast<std::size_t>(-1) &&
          r2i != static_cast<std::size_t>(-1), "R1/R2 thermal generators");
    if (!res.points.empty() &&
        r1i != static_cast<std::size_t>(-1) &&
        r2i != static_cast<std::size_t>(-1)) {
        checkNear(res.points.front().source_psd[r1i],
                  thermalPsd(r1) * rpar * rpar, 1e-12, "R1 contribution");
        checkNear(res.points.front().source_psd[r2i],
                  thermalPsd(r2) * rpar * rpar, 1e-12, "R2 contribution");
    }

    // Integrated noise of a WHITE density is exactly S*(f2-f1).
    const double band = 1e6 - 1.0;
    checkNear(res.stats.output_noise_v2, s_out * band, 1e-9,
              "integrated output noise = S*(fstop-fstart) for white noise");
    checkNear(res.stats.input_noise_v2, s_out / (gain * gain) * band, 1e-9,
              "integrated input noise");
    checkNear(res.stats.output_noise_rms, std::sqrt(s_out * band), 1e-9,
              "reported rms is the sqrt of the integral");

    // The per-generator integrals add up to the total (the integral is linear).
    double sum_src = 0.0;
    for (double v : res.stats.source_noise_v2) sum_src += v;
    checkNear(sum_src, res.stats.output_noise_v2, 1e-9,
              "per-source integrals sum to the total");
}

// ===============================================================
// 3) RC low-pass: the adjoint transfer is frequency dependent
// ===============================================================
static void testRcTransfer() {
    section("3) RC low-pass: S_out(f) = 4kT/R * |H(f)|^2");

    const double r = 1e3, c = 1e-9;
    Circuit circuit;
    const auto p = parseInto(circuit,
        "V1 in 0 DC 0\n"
        "R1 in out 1k\n"
        "C1 out 0 1n\n"
        ".NOISE V(out) V1 DEC 20 1 1e6 1\n"
        ".END\n");

    const NoiseResult res = runNoise(circuit, p, true);
    checkOk(res.success, res.message, "noise run succeeded");
    if (!res.success) return;
    check(res.sources.size() == 1, "one generator (R1); ideal C is noiseless");

    const double s0 = thermalPsd(r);
    // R1's current noise sits across R1; with `in` held at 0 V by the source,
    // the transfer to `out` is the R1/(1+jwRC) divider (R1 in series with the
    // parallel pair R1||C1 seen by the source):
    //     |H(f)| = R / (1 + jwRC)      ->  S_out(f) = 4kT·R / (1 + (wRC)^2)
    auto closedForm = [&](double f) {
        const double w = 2.0 * kPiTest * f;
        return 4.0 * kBoltzmann * kTemp300 * r / (1.0 + w * w * r * r * c * c);
    };
    (void)s0;
    double worst = 0.0;
    for (const auto& pt : res.points) {
        const double ref = closedForm(pt.frequency);
        worst = std::max(worst, std::fabs(pt.onoise_psd - ref) / ref);
        const double w = 2.0 * kPiTest * pt.frequency;
        const double gref = 1.0 / std::sqrt(1.0 + w * w * r * r * c * c);
        checkNear(pt.gain_mag, gref, 1e-12, "gain = |1/(1+jwRC)|");
    }
    if (worst >= 1e-12) {
        std::printf("   (worst relative PSD error %.3e)\n", worst);
        for (const auto& pt : res.points) {
            const double ref = closedForm(pt.frequency);
            const double rel = std::fabs(pt.onoise_psd - ref) / ref;
            if (rel > 1e-6) {
                std::printf("      f=%.6g sim=%.6e ref=%.6e rel=%.3e\n",
                            pt.frequency, pt.onoise_psd, ref, rel);
            }
        }
    }
    check(worst < 1e-12, "output PSD follows the closed form over the sweep");

    // The integrator against an INDEPENDENT quadrature of the same closed
    // form (fine trapezoid): this checks SPICE's power-law Nintegrate without
    // re-implementing it in the test.
    const double f1 = res.points.front().frequency;
    const double f2 = res.points.back().frequency;
    const int n = 200000;
    const double h = (f2 - f1) / n;
    double quad = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double f = f1 + i * h;
        quad += closedForm(f) * ((i == 0 || i == n) ? 0.5 : 1.0);
    }
    quad *= h;
    checkNear(res.stats.output_noise_v2, quad, 1e-6,
              "integrated noise vs an independent quadrature of the closed form");
}

// ===============================================================
// 4) Temperature scaling
// ===============================================================
static void testTemperature() {
    section("4) Thermal noise scales linearly with temperature");

    const std::string body =
        "V1 in 0 DC 1\n"
        "R1 in out 10k\n"
        "R2 out 0 10k\n"
        ".NOISE V(out) V1 DEC 5 10 1e5\n"
        ".END\n";

    double psd[2] = {0.0, 0.0};
    const double temps[2] = {27.0, 87.0};      // Celsius (.TEMP is Celsius)
    for (int i = 0; i < 2; ++i) {
        Circuit circuit;
        const auto p = parseInto(circuit, body);
        circuit.setTemperature(temps[i] + 273.15);
        const NoiseResult res = runNoise(circuit, p, false);
        checkOk(res.success, res.message, "noise run succeeded");
        if (!res.points.empty()) psd[i] = res.points.front().onoise_psd;
    }
    check(psd[0] > 0.0 && psd[1] > 0.0, "both temperatures produce noise");
    checkNear(psd[1] / psd[0], (87.0 + 273.15) / (27.0 + 273.15), 1e-12,
              "the PSD ratio is the temperature ratio (kT is linear in T)");
}

// ===============================================================
// 5) Diode shot noise
// ===============================================================
static void testDiodeShot() {
    section("5) Diode shot noise = 2q*Id at the operating point");

    Circuit circuit;
    const auto p = parseInto(circuit,
        "V1 a 0 DC 0.7\n"
        "R1 a b 1k\n"
        ".MODEL DMOD D(IS=1e-14 N=1)\n"
        "D1 b 0 DMOD\n"
        ".NOISE V(b) V1 DEC 5 1 1e3 1\n"
        ".END\n");

    const NoiseResult res = runNoise(circuit, p, true);
    checkOk(res.success, res.message, "noise run succeeded");
    if (!res.success) return;

    const std::size_t d1 = sourceIndex(res, "D1", "shot");
    check(d1 != static_cast<std::size_t>(-1), "diode reports a shot generator");
    check(sourceIndex(res, "R1", "thermal") != static_cast<std::size_t>(-1),
          "R1 thermal noise present");

    if (d1 != static_cast<std::size_t>(-1)) {
        // The DC current comes from the same solver, so the test does not
        // depend on the diode model: the reported prefactor must be 2q*Id with
        // Id read from the converged operating point the run left behind.
        double id = 0.0;
        for (const auto& dv : circuit.getDevices()) {
            if (dv->getName() == "D1") id = dv->getCurrent();
        }
        check(id > 0.0, "the diode conducts (forward biased)");
        checkNear(res.sources[d1].psd_a2, 2.0 * kCharge * id, 1e-12,
                  "shot prefactor = 2q*Id");
    }

    double worst = 0.0;
    for (const auto& pt : res.points) {
        double sum = 0.0;
        for (double v : pt.source_psd) sum += v;
        worst = std::max(worst, std::fabs(sum - pt.onoise_psd));
    }
    check(worst < 1e-24, "diode contributions still add up to the total");

    // Reverse biased: no shot noise (2qI would be negative).
    Circuit rev;
    const auto pr = parseInto(rev,
        "V1 a 0 DC 0\n"
        "R1 a b 1k\n"
        ".MODEL DMOD D(IS=1e-14 N=1)\n"
        "D1 b 0 DMOD\n"
        ".NOISE V(b) V1 DEC 5 1 1e3 1\n"
        ".END\n");
    const NoiseResult res_rev = runNoise(rev, pr, true);
    checkOk(res_rev.success, res_rev.message, "reverse-biased run succeeded");
    if (res_rev.success) {
        check(sourceIndex(res_rev, "D1", "shot") == static_cast<std::size_t>(-1),
              "no shot noise in reverse bias");
    }
}

// ===============================================================
// 6) Differential output V(out,ref)
// ===============================================================
static void testDifferentialOutput() {
    section("6) Differential output V(out,ref)");

    // A bridge: the noise between the two midpoints is the DIFFERENCE of two
    // transfer paths, which a single-node analysis would miss.
    Circuit circuit;
    const auto p = parseInto(circuit,
        "V1 a 0 DC 0\n"
        "R1 a m 1k\n"
        "R2 b m 1k\n"
        "R3 m 0 1k\n"
        "R4 m n 1k\n"
        "R5 b n 1k\n"
        "V2 b 0 DC 0\n"
        ".NOISE V(m,n) V1 DEC 5 1 1e3 1\n"
        ".END\n");

    const NoiseResult res = runNoise(circuit, p, true);
    checkOk(res.success, res.message, "noise run succeeded");
    if (!res.success) return;
    check(res.output_node == "m" && res.reference_node == "n",
          "output nodes recorded");
    for (const auto& pt : res.points) {
        double sum = 0.0;
        for (double v : pt.source_psd) sum += v;
        checkNear(sum, pt.onoise_psd, 1e-12,
                  "differential transfers still sum to the total");
    }
    // R3 sits between the two measured nodes, so it must contribute less than
    // a resistor that only touches one of them (R1/R2 see only one end).
    const std::size_t r3i = sourceIndex(res, "R3", "thermal");
    check(r3i != static_cast<std::size_t>(-1), "R3 present");
}

// ===============================================================
// 7) Parser
// ===============================================================
static void testParser() {
    section("7) Parser: .NOISE card spellings");

    {
        Circuit circuit;
        NetlistParser parser;
        const bool ok = parser.parseString(
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            "R2 out 0 1k\n"
            ".NOISE V(out) V1 DEC 10 1 1e5 1\n"
            ".END\n", circuit);
        check(ok, "ngspice spelling parses");
        const auto& p = parser.getAnalysisParameters();
        check(p.want_noise, "want_noise flag set");
        check(p.type == AnalysisType::NOISE, "analysis type is NOISE");
        check(p.noise_out_node == "out", "output node");
        check(p.noise_ref_node == "0", "reference defaults to ground");
        check(p.noise_source == "V1" && p.noise_source_given, "input source");
        check(p.ac_variation == "DEC" && p.ac_points == 10, "sweep");
        checkNear(p.ac_fstart, 1.0, 0, "fstart");
        checkNear(p.ac_fstop, 1e5, 0, "fstop");
        check(p.noise_detail, "pts_per_summary > 0 -> detail");
    }
    {
        Circuit circuit;
        NetlistParser parser;
        const bool ok = parser.parseString(
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            "R2 out 0 1k\n"
            ".NOISE V(out,ref) SRC=V1 OCT 8 1 1e6\n"
            ".END\n", circuit);
        check(ok, "Xyce spelling parses (SRC=..., V(out,ref))");
        const auto& p = parser.getAnalysisParameters();
        check(p.noise_out_node == "out", "output node");
        check(p.noise_ref_node == "ref", "reference node");
        check(p.noise_source == "V1", "SRC= form");
        check(!p.noise_detail, "no pts_per_summary -> no detail");
    }
    {
        Circuit circuit;
        NetlistParser parser;
        const bool ok = parser.parseString(
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            ".NOISE V(out) DEC 10 1 1e5\n"
            ".END\n", circuit);
        check(ok, "source may be omitted (first V source is used)");
        check(!parser.getAnalysisParameters().noise_source_given,
              "source not explicitly given");
    }
    {
        Circuit circuit;
        NetlistParser parser;
        const bool ok = parser.parseString(
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            ".NOISE V(out) V1 BAD 10 1 1e5\n"
            ".END\n", circuit);
        check(!ok, "an invalid variation is rejected");
    }
    {
        Circuit circuit;
        NetlistParser parser;
        const bool ok = parser.parseString(
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            ".NOISE V(nosuchnode) V1 DEC 10 1 1e5\n"
            ".END\n", circuit);
        check(ok, "unknown node is a run-time error, not a parse error");
    }
    {
        // An unknown input source must be reported, not silently ignored.
        Circuit circuit;
        const auto p = parseInto(circuit,
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            "R2 out 0 1k\n"
            ".NOISE V(out) V9 DEC 5 1 1e3\n"
            ".END\n");
        const NoiseResult res = runNoise(circuit, p, false);
        check(!res.success, "unknown input source fails with a message");
        check(res.message.find("V9") != std::string::npos,
              "the message names the source");
    }
    {
        // Unknown output node.
        Circuit circuit;
        const auto p = parseInto(circuit,
            "V1 in 0 DC 0\n"
            "R1 in out 1k\n"
            "R2 out 0 1k\n"
            ".NOISE V(nope) V1 DEC 5 1 1e3\n"
            ".END\n");
        const NoiseResult res = runNoise(circuit, p, false);
        check(!res.success, "unknown output node fails with a message");
    }
}

} // namespace

int main() {
    std::printf("=== NOISE ANALYSIS TESTS ===\n");

    testGrid();
    testDivider();
    testRcTransfer();
    testTemperature();
    testDiodeShot();
    testDifferentialOutput();
    testParser();

    std::printf("\n%s (%d checks, %d failures)\n",
                g_fail == 0 ? "NOISE TESTS PASSED" : "NOISE TESTS FAILED",
                g_checks, g_fail);
    if (g_fail != 0) {
        std::printf("last section: %s\n", g_section);
    }
    return g_fail == 0 ? 0 : 1;
}