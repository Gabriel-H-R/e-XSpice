/**
 * @author Gabriel Henrique Silva
 * @file ACAnalysis.h
 * @brief Small-signal AC analysis (.AC): one operating point, then a linear
 *        solve per frequency.
 *
 * Xyce-compatible subset (Reference Guide 7.8, section 2.2):
 *
 *   .AC <LIN|DEC|OCT> <points> <fstart> <fstop>      (frequency sweep)
 *
 * Sources carry an independent small-signal spec (`AC <mag> [phase]`), so
 * the same netlist can have a DC bias and an AC excitation.
 *
 * The formulation
 * ---------------
 * The simulator works with the DAE  C(x)·x' + F(x) = b(t).  Linearizing at
 * the DC operating point and going to the frequency domain gives
 *
 *     ( dF/dx + j·w·dQ/dx ) · X  =  B_ac(w)
 *
 * where G = dF/dx and C = dQ/dx are EXACTLY the two Jacobians the transient
 * linearization already builds at the operating point.  So the AC analysis
 * adds no device physics: it needs one new device hook (loadACExcitation,
 * the phasor on the right-hand side) and reuses the existing matrix
 * machinery for everything else.
 *
 * How this improves on Xyce's implementation
 * ------------------------------------------
 * Xyce (AnalysisPKG/N_ANP_AC.C) builds the AC system as a 2x2 REAL block
 * matrix
 *     [  G   -w·C ] [Vr]   [Br]
 *     [  w·C   G  ] [Vi] = [Bi]
 * i.e. a 2N x 2N real system: twice the factorization work (8/3·N³ real flops
 * instead of 4/3·N³ for a complex LU) and twice the memory. This engine
 * factorizes the complex matrix directly (ComplexSparseLU), which:
 *
 *  - reuses the DC sparsity pattern and computes the symbolic factor ONCE;
 *    only the numeric values are refilled per frequency (the device
 *    linearization is evaluated a single time for the whole sweep),
 *  - stores values as two real arrays (SoA) — better vectorization, no
 *    std::complex temporaries,
 *  - runs FREQUENCIES in parallel (OpenMP over the sweep, optionally
 *    distributed across MPI ranks with setPointDistribution), instead of
 *    one factorization at a time in a serial loop. Frequencies are
 *    independent, and each factorization is self-contained, so the result
 *    does not depend on how the work was split: a parallel run is
 *    bit-for-bit identical to the serial one (unlike the DC sweep, whose
 *    continuation chain makes it order dependent).
 *
 * Deliberate limitation: the factorization keeps the natural ordering and
 * does no pivoting (that is what makes the static pattern and the reuse
 * possible), with a relative pivot threshold plus diagonal perturbation as
 * the safeguard. Circuit matrices (conductances of a passive network plus
 * linearized devices) are diagonally dominant, so this is the right trade;
 * a pathological matrix is reported instead of silently mis-solved.
 */

#ifndef AC_ANALYSIS_H
#define AC_ANALYSIS_H

#include "../core/Circuit.h"
#include "../solver/LinearSolver.h"
#include "DaeSystem.h"

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief Frequency sweep spacing (same spelling as the .DC sweep).
 */
enum class ACVariation {
    LIN,   // .AC LIN <points> <fstart> <fstop>   (points in total)
    DEC,   // .AC DEC <points> ...                 (points per decade)
    OCT    // .AC OCT <points> ...                 (points per octave)
};

/**
 * @brief One frequency point of the result.
 *
 * The unknown vector X holds node voltages in the first numNodeRows() entries
 * and branch currents (V sources, inductors) in the rest, so the phasor of
 * I(V1) is simply the entry of the corresponding extra row — no separate
 * device query is needed.
 */
struct ACDataPoint {
    double frequency = 0.0;
    std::vector<double> mag;       // |X[row]|            per node row
    std::vector<double> phase_deg;  // arg(X[row]) degrees
    std::vector<double> db;         // 20·log10(|X[row]|)  per node row
    std::vector<double> branch_mag;      // |I| per extra unknown (source/inductor)
    std::vector<double> branch_phase_deg;
    bool valid = true;             // false when this frequency failed
};

/**
 * @brief AC sweep accounting (performance transparency).
 */
struct ACStats {
    long long frequencies = 0;        // points solved
    int  factorizations = 0;          // complex LU factorizations
    long long factor_flops = 0;       // accumulation flops of those
    int  newton_iterations = 0;       // operating-point solve
    int  failed_frequencies = 0;
    bool pivot_perturbed = false;     // safeguard triggered (should be rare)
    int  ranks_used = 1;              // 1 = local parallel, >1 = split
    double elapsed_seconds = 0.0;
};

/**
 * @brief Result of the AC analysis.
 */
struct ACResult {
    bool success = false;
    std::string message;
    std::string variation = "DEC";
    std::vector<ACDataPoint> points;   // decimated sample (setMaxRetainedPoints)
    ACStats stats;
};

/**
 * @brief AC analysis driver.
 *
 * Lifecycle: construct -> configure -> run(). The operating point is solved
 * once; the (G, C) snapshot is taken once; every frequency is then an
 * independent complex solve.
 */
class ACAnalysis {
public:
    explicit ACAnalysis(Circuit& circuit);
    ~ACAnalysis();

    // ---- sweep definition ----
    void setVariation(ACVariation v) { variation_ = v; }
    void setPoints(int n) { points_ = n; }
    void setFrequencyRange(double fstart, double fstop) {
        fstart_ = fstart; fstop_ = fstop;
    }

    // ---- solver knobs ----
    void setMaxIterations(int n) { max_iterations_ = n; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool v) { verbose_ = v; }

    /**
     * @brief Splits the FREQUENCY sweep across MPI ranks
     *        (`.OPTIONS ACFREQPAR`).
     *
     * Unlike the DC sweep, this is not a continuation chain: frequencies are
     * independent and each rank factorizes the same replicated (G, C) with
     * the same code, so the answer is bit-identical to a serial run — the
     * only cost is the exchange of the resulting phasors (batched, one
     * gather per kExchangeChunk frequencies, not one per frequency).
     *
     * It pays off when a single factorization is significant (many unknowns
     * and/or many frequencies). For a tiny circuit the exchange dominates
     * and the default is faster.
     */
    void setPointDistribution(bool on) { point_distribution_ = on; }
    bool pointDistribution() const { return point_distribution_; }

    /**
     * @brief Upper bound on points kept in ACResult::points; the CSV stream
     *        (when set) still receives EVERY frequency. 0 = keep all.
     */
    void setMaxRetainedPoints(size_t n) { max_retained_ = n; }

    /**
     * @brief Streams every solved frequency as CSV rows while the sweep runs.
     * @param labels column labels for the node phasors (node order)
     */
    void setCsvStream(std::ostream* os, std::vector<std::string> labels);

    void setProgressCallback(std::function<void(size_t, size_t, double)> cb) {
        progress_ = std::move(cb);
    }

    /**
     * @brief Builds the frequency grid (pure; exposed for tests).
     *
     * LIN spreads `points` values linearly from fstart to fstop; DEC/OCT use
     * `points` values per decade/octave and always INCLUDE both endpoints,
     * exactly like the .DC sweep and Xyce.
     */
    static std::vector<double> buildGrid(ACVariation variation, int points,
                                         double fstart, double fstop);

    /**
     * @brief Runs: operating point -> linearization -> frequency sweep.
     * @param use_mpi enables MPI distribution (frequencies) when set
     */
    ACResult run(bool use_mpi = true);

private:
    Circuit& circuit_;
    ACVariation variation_ = ACVariation::DEC;
    int    points_ = 10;
    double fstart_ = 1.0;
    double fstop_  = 1e6;

    int    max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool   verbose_ = false;
    bool   point_distribution_ = false;
    size_t max_retained_ = 0;

    std::ostream* csv_ = nullptr;
    std::vector<std::string> csv_labels_;
    bool csv_header_written_ = false;
    std::streamsize csv_saved_precision_ = 0;

    std::function<void(size_t, size_t, double)> progress_;

    static constexpr int kMaxFrequencies_ = 1000000;
    static constexpr int kExchangeChunk_ = 256;  // frequencies per gather
};

} // namespace e_XSpice

#endif // AC_ANALYSIS_H