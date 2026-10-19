/**
 * @author Gabriel Henrique Silva
 * @file NoiseAnalysis.h
 * @brief Noise analysis (.NOISE): output/input noise spectral densities and
 *        the integrated noise, via ONE adjoint solve per frequency.
 *
 * Card (SPICE3 / ngspice / Xyce compatible)
 * ---------------------------------------
 *     .NOISE V(out[,ref]) <src> (DEC|OCT|LIN) <pts> <fstart> <fstop>
 *           [pts_per_summary]
 *
 *   - `V(out[,ref])` is the output: the noise of the voltage between the two
 *     nodes (ref defaults to ground).
 *   - `src` names the independent source the input noise is referred to. It may
 *     also be written `SRC=V1` (Xyce's spelling). If omitted, the first
 *     independent voltage source of the circuit is used.
 *   - `pts/fstart/fstop` are the .AC sweep parameters.
 *   - `pts_per_summary` (optional, > 0) asks for the per-generator
 *     contribution breakdown.
 *
 * Method
 * ------
 * At the DC operating point the small-signal system is A = G + jwC. Every noise
 * generator of the circuit is a Norton current pair with a known PSD; the noise
 * the output sees is
 *
 *      S_out(w) = sum_i |H_i(w)|^2 · S_i(w),      H_i = A^-1·(e_pos - e_neg)
 *
 * Solving that per generator would cost ONE LINEAR SOLVE PER GENERATOR per
 * frequency — the classic O(n_sources x N) blow-up. By symmetry of A:
 *
 *      A^T · x* = e_out     =>     H_i = x*[pos] - x*[neg]
 *
 * so a single adjoint solve yields the transfer of ALL generators at once
 * (this is what Xyce does too — see AnalysisPKG/N_ANP_NOISE.C,
 * solveAdjointNOISE_ — and what SPICE3 did originally). On top of that this
 * engine adds:
 *
 *  - ONE solve instead of two: the AC gain V(out)/V(src) needed for the
 *    input-referred noise comes from the SAME adjoint vector by reciprocity
 *    (the RHS of a voltage source is its branch unknown, so gain = x*[branch]),
 *    which also makes the noise analysis independent of the AC magnitudes in
 *    the netlist. Xyce needs an extra AC solve per frequency for the gain.
 *  - The frequency-independent part of every generator (4kT/R, 2qI, KF·I^AF…)
 *    is computed ONCE at the operating point; the sweep only applies the
 *    f^(-NF) factor. Xyce calls a virtual getNoiseSources() on every device at
 *    every frequency.
 *  - Frequencies run in parallel (OpenMP, optionally split across MPI ranks
 *    with ACFREQPAR), and results are bit-identical to the serial run because
 *    frequencies are independent and the matrix is replicated.
 *  - Output is streamed; only a decimated sample plus O(n_sources) integral
 *    accumulators stay in memory, instead of a full per-source spectrum table.
 *  - The complex factorization is the same ComplexSparseLU of the AC analysis
 *    (a native complex LU, half the flops and the memory of the 2N x 2N real
 *    block system Xyce builds for AC/noise), with the symbolic factorization
 *    done once for the whole sweep.
 *
 * Integrals follow SPICE3's Nintegrate (a piecewise power law in log-log
 * space), which is also Xyce's noiseIntegral(): the density between two sweep
 * points is fitted as a·f^k and integrated in closed form. Integrated values
 * are reported as V^2 (the integral of a PSD) and as V rms (what SPICE calls
 * "total noise").
 */

#ifndef NOISE_ANALYSIS_H
#define NOISE_ANALYSIS_H

#include "../core/Circuit.h"
#include "ACAnalysis.h"          // ACVariation (same sweep semantics)

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief One noise generator, with its terminals already mapped onto the
 *        ground-eliminated frequency system.
 */
struct NoiseSourceInfo {
    std::string device;
    std::string label;
    int  pos = -1;             // frequency-system row (or -1 = ground)
    int  neg = -1;
    double psd_a2 = 0.0;       // S(1 Hz) in A^2/Hz
    double f_exponent = 0.0;   // S(f) = psd_a2 · f^(-f_exponent)
    std::string column() const { return device + "_" + label + "_psd"; }
};

/**
 * @brief One frequency point of the result. Densities are POWER spectral
 *        densities (V^2/Hz); multiply by dB conversion with 10·log10.
 */
struct NoiseDataPoint {
    double frequency = 0.0;
    double onoise_psd = 0.0;      // V^2/Hz at the output
    double inoise_psd = 0.0;      // V^2/Hz referred to the input
    double gain_mag = 0.0;        // |V(out)/V(src)|
    double gain_phase_deg = 0.0;
    double onoise_amp = 0.0;      // V/sqrt(Hz) (what a designer reads)
    double inoise_amp = 0.0;      // V/sqrt(Hz)
    std::vector<double> source_psd;   // per generator, V^2/Hz at the output
    bool valid = true;
};

/**
 * @brief Noise sweep accounting + the integrated summary.
 */
struct NoiseStats {
    long long frequencies = 0;
    int  factorizations = 0;
    long long factor_flops = 0;
    int  newton_iterations = 0;
    int  failed_frequencies = 0;
    int  sources = 0;
    bool pivot_perturbed = false;
    int  ranks_used = 1;
    double elapsed_seconds = 0.0;

    double output_noise_v2 = 0.0;      // integral of the output PSD (V^2)
    double input_noise_v2 = 0.0;       // integral of the input PSD  (V^2)
    double output_noise_rms = 0.0;     // V
    double input_noise_rms = 0.0;      // V
    double peak_output_psd = 0.0;      // max of onoise over the sweep
    double peak_output_freq = 0.0;

    /** Per-generator integrated output noise (V^2), aligned with the source
     *  list of the result. */
    std::vector<double> source_noise_v2;
};

/**
 * @brief Result of the noise analysis.
 */
struct NoiseResult {
    bool success = false;
    std::string message;
    std::string variation = "DEC";
    std::string output_node;
    std::string reference_node = "0";
    std::string input_source;
    std::vector<NoiseSourceInfo> sources;
    std::vector<NoiseDataPoint> points;    // decimated sample
    NoiseStats stats;
};

class NoiseAnalysis {
public:
    explicit NoiseAnalysis(Circuit& circuit);
    ~NoiseAnalysis();

    // ---- sweep definition ----
    void setVariation(ACVariation v) { variation_ = v; }
    void setPoints(int n) { points_ = n; }
    void setFrequencyRange(double fstart, double fstop) {
        fstart_ = fstart; fstop_ = fstop;
    }

    // ---- output specification ----
    /** Output is V(out) (reference = ground). */
    void setOutputNode(const std::string& node);
    /** Output is V(out) - V(ref). */
    void setOutputNode(const std::string& node, const std::string& ref);

    /** Independent source the input noise is referred to (e.g. "V1").
     *  Empty = the first independent voltage source of the circuit. */
    void setInputSource(const std::string& name) { input_source_ = name; }

    /** Per-generator contribution breakdown (SPICE's pts_per_summary). */
    void setSourceDetail(bool on) { source_detail_ = on; }

    // ---- solver knobs ----
    void setMaxIterations(int n) { max_iterations_ = n; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool v) { verbose_ = v; }
    void setPointDistribution(bool on) { point_distribution_ = on; }
    bool pointDistribution() const { return point_distribution_; }

    /** Upper bound on points kept in NoiseResult::points; the CSV stream
     *  (when set) still receives EVERY frequency. 0 = keep all. */
    void setMaxRetainedPoints(size_t n) { max_retained_ = n; }

    void setCsvStream(std::ostream* os) { csv_ = os; }

    void setProgressCallback(std::function<void(size_t, size_t, double)> cb) {
        progress_ = std::move(cb);
    }

    /** Builds the frequency grid (same rule as the AC sweep). */
    static std::vector<double> buildGrid(ACVariation variation, int points,
                                         double fstart, double fstop);

    /**
     * @brief Runs: operating point -> generator list -> adjoint sweep.
     * @param use_mpi enables MPI frequency distribution when set
     */
    NoiseResult run(bool use_mpi = true);

private:
    Circuit& circuit_;
    ACVariation variation_ = ACVariation::DEC;
    int    points_ = 10;
    double fstart_ = 1.0;
    double fstop_  = 1e6;

    std::string output_node_;
    std::string reference_node_ = "0";
    std::string input_source_;
    bool source_detail_ = false;

    int    max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool   verbose_ = false;
    bool   point_distribution_ = false;
    size_t max_retained_ = 0;

    std::ostream* csv_ = nullptr;
    bool csv_header_written_ = false;
    std::streamsize csv_saved_precision_ = 0;

    std::function<void(size_t, size_t, double)> progress_;

    static constexpr int kMaxFrequencies_ = 1000000;
    static constexpr int kExchangeChunk_ = 256;
};

} // namespace e_XSpice

#endif // NOISE_ANALYSIS_H