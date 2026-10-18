/**
 * @author Gabriel Henrique Silva
 * @file DCSweep.h
 * @brief DC sweep analysis (.DC): one DAE setup, continuation chain of OPs.
 *
 * Xyce-compatible subset (Reference Guide 7.8, section 2.1.3):
 *
 *   .DC [LIN] <var> <start> <stop> <step>
 *   .DC DEC   <var> <start> <stop> <points-per-decade>
 *   .DC OCT   <var> <start> <stop> <points-per-octave>
 *   .DC <spec1> [[LIN|DEC|OCT] <var2> <start> <stop> <step>]   (nested)
 *
 * Sweepable variables: independent V and I sources, and resistors
 * (`.DC R1 0 3.5 0.05` in the Xyce guide). Each level of the nested pair is
 * one of those; the second specification is the INNER loop and restarts from
 * its own start value at every point of the outer loop, exactly like Xyce.
 *
 * Why this is not "one DCAnalysis per point"
 * -------------------------------------------
 *  - ONE DaeSystem::setup() for the whole sweep. Across sweep points only
 *    source values move: the sparsity pattern, the CSR slot assignment, the
 *    StateStore layout and the MPI device/row plan are invariant, so they
 *    are built once instead of N times. For a 1000-point sweep on a large
 *    circuit this is where the win is, not in the Newton loop.
 *  - Honest continuation: each point seeds Newton from the previous
 *    converged solution, so the walk follows the solution branch instead of
 *    re-discovering it from zero at every point.
 *  - Linear predictor: when the step ratio is regular, the seed is
 *    extrapolated from the last two accepted points (Xyce's continuation
 *    predictor), typically saving one Newton iteration per point.
 *  - Source-step refinement on failure (Xyce retries a failed bias point
 *    with source stepping instead of aborting): the failing step is bisected
 *    and the intermediate points are solved, which recovers sweeps that a
 *    single shot would lose.
 *  - Bounded memory: CSV rows are streamed as points are produced and the
 *    in-memory sample is decimated, so RAM is O(nodes) and not
 *    O(points x nodes).
 *  - Temperature as a sweep variable (`.DC TEMP <start> <stop> <step>` in
 *    Celsius, Xyce): one Circuit::setTemperature() per point re-scales the
 *    shared cards in place, so the SAME DaeSystem keeps solving — a
 *    temperature sweep reuses the setup, pattern and slots.
 *  - MPI: by default the sweep is deterministic and identical on every rank
 *    (no new collectives); device-parallel assembly and the row-distributed
 *    Krylov solve live inside DaeSystem. `.OPTIONS SWEEPPAR` splits the sweep
 *    POINTS across ranks instead — see setPointDistribution() for the exact
 *    trade-offs.
 */

#ifndef DC_SWEEP_H
#define DC_SWEEP_H

#include "../core/Circuit.h"
#include "../solver/LinearSolver.h"
#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief Sweep grid spacing.
 */
enum class DCSweepVariation {
    LIN,   // .DC [LIN] src start stop step
    DEC,   // .DC DEC src start stop points-per-decade
    OCT    // .DC OCT src start stop points-per-octave
};

/**
 * @brief One sweep level (variable + range).
 */
struct DCSweepSpec {
    std::string source;
    DCSweepVariation variation = DCSweepVariation::LIN;
    double start = 0.0;
    double stop  = 0.0;
    double step  = 0.0;   // LIN: step; DEC/OCT: points per decade/octave
    bool valid() const { return !source.empty(); }
};

/**
 * @brief One requested (retained) sweep point.
 */
struct DCSweepPoint {
    double sweep_value = 0.0;             // inner (or only) sweep value
    double outer_value = 0.0;             // outer sweep value (0 if not nested)
    std::vector<double> node_voltages;    // per Circuit::getNodes() order
    std::vector<double> branch_currents;  // per Circuit::getDevices() order
    int  newton_iterations = 0;
    bool converged = true;                // false: NOFAIL row (holds last)
};

/**
 * @brief Sweep accounting (performance transparency + NOFAIL report).
 */
struct DCSweepStats {
    long long points_solved = 0;      // requested points attempted
    long long newton_iterations = 0;  // summed over every solve
    long long linear_solves = 0;      // J*dx systems
    long long predictor_hits = 0;     // points seeded by extrapolation
    int  retry_solves = 0;            // source-step refinement solves
    int  failed_points = 0;           // points that never converged
    int  max_point_iterations = 0;
    int  ranks_used = 1;              // 1 = sequential chain, >1 = points split
    double elapsed_seconds = 0.0;
};

/**
 * @brief DC sweep result.
 */
struct DCSweepResult {
    bool success = false;
    std::string message;
    std::string sweep_source;             // outer (or only) variable name
    std::string sweep_source2;            // inner variable (empty if flat)
    std::string variation  = "LIN";
    std::string variation2 = "LIN";
    bool nested = false;
    std::vector<DCSweepPoint> points;     // decimated sample (see setMaxRetainedPoints)
    DCSweepStats stats;
};

/**
 * @brief A quantity a .DC sweep can drive: independent V/I source, resistor,
 *        or the circuit temperature (Xyce `.DC ... TEMP ...`).
 *
 * Adapters over the devices and over Circuit::setTemperature(); the sweep
 * loop only sees get()/set(), so no device type ever leaks into the analysis
 * code. TEMP values are degrees Celsius (SPICE/Xyce convention); the adapter
 * converts to Kelvin for the circuit.
 */
class SweepTarget {
public:
    virtual ~SweepTarget() = default;
    virtual double      get() const = 0;
    virtual void        set(double value) = 0;
    virtual const char* kind() const = 0;
};

/**
 * @brief DC sweep driver on top of ONE shared DaeSystem.
 */
class DCSweep {
public:
    explicit DCSweep(Circuit& circuit);
    ~DCSweep();

    // ---- sweep definition ----
    void setSweep(const DCSweepSpec& spec) { spec_ = spec; }
    void setInnerSweep(const DCSweepSpec& spec) { spec2_ = spec; }
    void setSource(const std::string& name) { spec_.source = name; }
    void setVariation(DCSweepVariation v) { spec_.variation = v; }
    void setRange(double start, double stop, double step_or_points) {
        spec_.start = start; spec_.stop = stop; spec_.step = step_or_points;
    }

    // ---- solver knobs ----
    void setSolverType(SolverType type) { solver_type_ = type; }
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool verbose) { verbose_ = verbose; }

    /**
     * @brief Keep sweeping after a point fails (`.OPTIONS NOFAIL`).
     *        The failed row is emitted flagged with the last converged
     *        solution. Default false: a failed point aborts the sweep.
     */
    void setContinueOnFail(bool on) { continue_on_fail_ = on; }

    /**
     * @brief Source-step refinement depth for a failed point: up to N
     *        bisections of the failing step (each one is a real solve).
     */
    void setMaxRetries(int n) { max_retries_ = std::max(0, n); }

    /**
     * @brief Linear predictor for the Newton seed (default on).
     */
    void setPredictor(bool on) { predictor_ = on; }

    /**
     * @brief Distributes the sweep POINTS across MPI ranks
     *        (`.OPTIONS SWEEPPAR`).
     *
     * Why this is opt-in and different from device distribution: a .DC
     * sweep is a dependent chain (each point warm-starts from the previous
     * one and from the device limiter history), so splitting it means every
     * rank runs its OWN chain over its own subset of the grid, replicating
     * the circuit. Consequences, stated plainly:
     *  - the default (off) keeps the sequential chain: identical on every
     *    rank, bit-for-bit reproducible CSV, no extra communication, and
     *    parallelism comes from DaeSystem (device assembly + distributed
     *    Krylov).
     *  - with it on, ranks converge from different seeds, so the answer is
     *    the same root but not bit-identical to the sequential run. The
     *    deviation is the path-dependent numerical noise of the nonlinear
     *    chain itself, bounded by the Newton criterion (RELTOL): measured
     *    ~8e-6 V on a diode clipper and ~1.6e-3 V on a 400-section diode
     *    ladder — the same order as the difference between two SEQUENTIAL
     *    runs of that ladder at different grid steps (9e-4 V).
     *    Results travel in chunks (one gather per 256 points, not per
     *    point) and the master writes the CSV in sweep order; points and
     *    CSV live on the master, while the outcome (success/message) and
     *    the statistics are broadcast so every rank reports the same thing.
     *  - it pays off when the PER-POINT SOLVE is expensive (many devices,
     *    many points): measured 1.5x on 4 ranks for a 400-node diode ladder
     *    with 101 sweep points. On a 3-device circuit the batched gather and
     *    the MPI startup dominate and the mode is SLOWER than the
     *    sequential chain — leave it off there, where the default also
     *    stays bit-reproducible.
     *
     * Nested sweeps distribute the OUTER loop (each outer row is an
     * independent chain of inner points).
     */
    void setPointDistribution(bool on) { point_distribution_ = on; }
    bool pointDistribution() const { return point_distribution_; }

    /**
     * @brief Upper bound on points kept in DCSweepResult::points. The CSV
     *        stream (when set) still receives EVERY point, so this only caps
     *        RAM: 0 = keep everything (tests / small sweeps).
     */
    void setMaxRetainedPoints(size_t n) { max_retained_ = n; }

    /**
     * @brief Streams every solved point as a CSV row while the sweep runs.
     * @param labels column headers for the node voltages, in
     *        Circuit::getNodes() order (as produced by main's nodeLabels).
     */
    void setCsvStream(std::ostream* os, std::vector<std::string> labels);

    /**
     * @brief Progress hook, called with (points_done, points_total, value) —
     *        the transient engine has the equivalent one.
     */
    void setProgressCallback(std::function<void(size_t, size_t, double)> cb) {
        progress_ = std::move(cb);
    }

    /**
     * @brief Builds a sweep grid for one level (pure; exposed for tests).
     */
    std::vector<double> buildGrid(const DCSweepSpec& spec) const;
    std::vector<double> buildGrid() const { return buildGrid(spec_); }

    /**
     * @brief Resolves a netlist name to a sweepable target (V/I source or R).
     * @param error filled with a user-facing message when unsupported.
     */
    static std::unique_ptr<SweepTarget> makeTarget(Circuit& circuit,
                                                   const std::string& name,
                                                   std::string& error);

    /**
     * @brief Runs the sweep.
     * @param use_mpi enables MPI distribution inside the shared DaeSystem
     *        (and point distribution when setPointDistribution(true))
     */
    DCSweepResult run(bool use_mpi = true);

private:
    Circuit& circuit_;
    DCSweepSpec spec_;    // outer (or only) sweep level
    DCSweepSpec spec2_;   // inner sweep level (invalid = flat sweep)

    SolverType solver_type_ = SolverType::BICGSTAB;
    int  max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool verbose_ = false;

    bool continue_on_fail_ = false;
    int  max_retries_ = 3;
    bool predictor_ = true;
    bool point_distribution_ = false;
    size_t max_retained_ = 0;

    std::ostream* csv_ = nullptr;
    std::vector<std::string> csv_labels_;
    bool csv_header_written_ = false;
    std::streamsize csv_saved_precision_ = 0;

    std::function<void(size_t, size_t, double)> progress_;

    static constexpr int kMaxPoints_ = 100000;   // grid hard cap per level
    // Predictor guards: only extrapolate when the previous point was cheap
    // (a hard point has a jagged path) and the step ratio is regular.
    static constexpr int kPredictorMaxIters_ = 6;
    static constexpr double kMaxPredictorRatio_ = 2.0;
};

} // namespace e_XSpice

#endif // DC_SWEEP_H