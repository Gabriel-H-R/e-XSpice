/**
 * @author Gabriel Henrique Silva
 * @file DCSweep.cpp
 * @brief DC sweep: one DaeSystem setup, continuation chain, streamed output.
 *
 * Traversal is organized around CHAIN UNITS:
 *   - flat sweep   : one unit per point
 *   - nested sweep : one unit per OUTER value (its inner sweep is a chain)
 * A unit is solved by SweepEngine, which owns the continuation state, the
 * predictor and the failure ladder. The two drivers differ only in who owns
 * which units: runSequential() walks every unit on every rank (identical,
 * bit-reproducible), runDistributed() splits the units across ranks and
 * exchanges the resulting points per point.
 */

#include "DCSweep.h"
#include "DaeSystem.h"
#include "../devices/VoltageSource.h"
#include "../devices/CurrentSource.h"
#include "../devices/Resistor.h"
#include "../devices/SourceWaveform.h"
#include "../parallel/MPIManager.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

namespace e_XSpice {

namespace {

// ---- sweep target adapters -------------------------------------------------
// One tiny class per supported variable kind; the sweep loop itself never
// mentions a device type (same discipline as the DAE device contract).

class VSourceTarget : public SweepTarget {
public:
    explicit VSourceTarget(VoltageSource* d) : d_(d) {}
    double get() const override { return d_->waveform().value(0.0); }
    void set(double v) override { d_->setWaveform(SourceWaveform::dc(v)); }
    const char* kind() const override { return "voltage source"; }
private:
    VoltageSource* d_;   // non-owning (devices are RCP-owned by Circuit)
};

class ISourceTarget : public SweepTarget {
public:
    explicit ISourceTarget(CurrentSource* d) : d_(d) {}
    double get() const override { return d_->waveform().value(0.0); }
    void set(double v) override { d_->setWaveform(SourceWaveform::dc(v)); }
    const char* kind() const override { return "current source"; }
private:
    CurrentSource* d_;
};

class ResistorTarget : public SweepTarget {
public:
    explicit ResistorTarget(Resistor* d) : d_(d) {}
    double get() const override { return d_->getEffectiveResistance(); }
    void set(double v) override { d_->setEffectiveResistance(v); }
    const char* kind() const override { return "resistor"; }
private:
    Resistor* d_;
};

/**
 * @brief `.DC ... TEMP ...`: the circuit temperature, in Celsius.
 *
 * Circuit::setTemperature() re-scales every model card that did not pin
 * its own TEMP and refreshes the devices that precompute temperature
 * scalars, so the SAME DaeSystem (pattern, slots, states) keeps solving —
 * a temperature sweep costs one re-scaling per point, not a new setup.
 */
class TemperatureTarget : public SweepTarget {
public:
    explicit TemperatureTarget(Circuit* c) : c_(c) {}
    double get() const override { return c_->temperature() - Const::CtoK; }
    void set(double v) override { c_->setTemperature(v + Const::CtoK); }
    const char* kind() const override { return "circuit temperature"; }
private:
    Circuit* c_;
};

const char* variationName(DCSweepVariation v) {
    switch (v) {
        case DCSweepVariation::DEC: return "DEC";
        case DCSweepVariation::OCT: return "OCT";
        default: return "LIN";
    }
}

// Restores the swept quantities when the sweep leaves (error path or end).
struct ScopeGuard {
    std::function<void()> fn;
    ~ScopeGuard() { if (fn) fn(); }
};

/**
 * @brief Solver-side configuration of a chain walk (plain data, so the
 *        engine stays independent from the DCSweep option setters).
 */
struct EngineConfig {
    bool predictor = true;
    bool continue_on_fail = false;
    bool verbose = false;
    int  max_retries = 3;
    // A point that the predictor may extrapolate for: only when the previous
    // point was cheap and the step ratio is regular.
    int  predictor_max_iters = 6;
    double predictor_max_ratio = 2.0;
};

/**
 * @brief One solved point, ready to be reported or exchanged.
 */
struct PointRow {
    size_t index = 0;              // global position in the sweep
    double value = 0.0;            // inner (or only) sweep value
    double outer = 0.0;            // outer sweep value
    int    iters = 0;
    bool   converged = true;
    std::vector<double> nodes;     // node rows of the solution
    std::vector<double> currents;  // device currents (retained points only)
};

/**
 * @brief Chain walk over units: continuation, predictor, failure ladder.
 *
 * Owns everything that makes a point cheap (warm start) or recoverable
 * (source stepping); the drivers only decide which units to walk.
 */
class SweepEngine {
public:
    SweepEngine(DaeSystem& system, Circuit& circuit, SweepTarget& primary,
                SweepTarget* secondary, const EngineConfig& cfg)
        : system_(system)
        , circuit_(circuit)
        , primary_(primary)
        , secondary_(secondary)
        , cfg_(cfg)
        , row_(static_cast<size_t>(std::max(0, system.numNodeRows())))
    {}

    /**
     * @brief Solves one chain unit.
     *
     * @param outer_grid  outer values (flat sweep: the point values)
     * @param inner_grid  inner values (empty = flat sweep)
     * @param unit        unit index (row of outer_grid, or point index)
     * @param base_index  global index of the unit's first point
     * @param out         solved points, appended in sweep order
     * @return false when the unit aborted the sweep (failure, no NOFAIL)
     */
    bool solveUnit(const std::vector<double>& outer_grid,
                   const std::vector<double>& inner_grid,
                   size_t unit, size_t base_index,
                   std::vector<PointRow>& out) {
        if (aborted_) return false;
        const bool nested = !inner_grid.empty();
        const double outer_value = outer_grid[unit];
        if (secondary_) secondary_->set(outer_value);

        const size_t n_in = nested ? inner_grid.size() : 1;
        for (size_t ii = 0; ii < n_in; ++ii) {
            const double v = nested ? inner_grid[ii] : outer_grid[unit];
            const size_t index = base_index + ii;

            // ---- seed: continuation + linear extrapolation ----------------
            seed_.clear();
            if (!x_last_.empty()) {
                seed_ = x_last_;
                if (cfg_.predictor && !x_ante_.empty()
                    && last_iters_ <= cfg_.predictor_max_iters) {
                    const double dv_prev = last_v_ - ante_v_;
                    const double dv_next = v - last_v_;
                    if (std::abs(dv_prev) > 0.0) {
                        const double ratio = dv_next / dv_prev;
                        // A wild ratio means the solution path is not locally
                        // linear (decade jump, restart of the inner loop):
                        // extrapolating there is worse than plain continuation.
                        if (std::isfinite(ratio)
                            && std::abs(ratio) <= cfg_.predictor_max_ratio) {
                            for (size_t k = 0; k < seed_.size(); ++k) {
                                seed_[k] = x_last_[k]
                                         + (x_last_[k] - x_ante_[k]) * ratio;
                            }
                            ++predictor_hits_;
                        }
                    }
                }
            }

            NLSolveResult nl = solveAt(v, seed_);
            bool ok = nl.converged;

            // ---- failure ladder (Xyce retries a failed bias point with
            //      source stepping rather than aborting the sweep):
            //      1. cold start, when there is no chain to continue from;
            //      2. source-step refinement: bisect the failing step, each
            //         intermediate point being a real operating point;
            //      3. the target again, warm-started from the refined point.
            if (!ok && !chain_) {
                ++retry_solves_;
                cold_.assign(static_cast<size_t>(system_.size()), 0.0);
                system_.seedInitialGuess(cold_);
                nl = solveAt(v, cold_);
                ok = nl.converged;
            }
            if (!ok && chain_ && cfg_.max_retries > 0) {
                double lo = last_v_;
                double hi = v;
                for (int k = 0; k < cfg_.max_retries; ++k) {
                    const double mid = 0.5 * (lo + hi);
                    ++retry_solves_;
                    if (solveAt(mid, x_last_).converged) lo = mid;
                    else                                  hi = mid;
                }
                ++retry_solves_;
                nl = solveAt(v, x_last_);
                ok = nl.converged;
            }

            ++points_solved_;
            if (!ok) {
                ++failed_points_;
                if (!cfg_.continue_on_fail) {
                    aborted_ = true;
                    abort_value_ = v;
                    abort_outer_ = outer_value;
                    abort_index_ = index;
                    abort_message_ = nl.message;
                    return false;
                }
            }

            // ---- row: node rows of the solution vector (extractSolution
            //      maps them 1:1 onto Circuit::getNodes(), so this is a
            //      straight copy with no virtual dispatch per node).
            const std::vector<double>& src =
                ok ? system_.solution() : (chain_ ? x_last_ : cold_);
            if (ok || chain_) {
                std::copy_n(src.begin(), row_.size(), row_.begin());
            } else {
                std::fill(row_.begin(), row_.end(), 0.0);
            }

            PointRow r;
            r.index = index;
            r.value = v;
            r.outer = outer_value;
            r.iters = last_iters_;
            r.converged = ok;
            r.nodes = row_;
            out.push_back(std::move(r));
        }
        return true;
    }

    /**
     * @brief Fills the device currents of a reported point (retained points
     *        only: the exchange payload stays small in distributed mode).
     */
    void addCurrents(PointRow& r) const {
        r.currents.reserve(circuit_.getDevices().size());
        for (const auto& d : circuit_.getDevices()) {
            r.currents.push_back(d->getCurrent());
        }
    }

    bool aborted() const { return aborted_; }
    double abortValue() const { return abort_value_; }
    double abortOuter() const { return abort_outer_; }
    size_t abortIndex() const { return abort_index_; }
    const std::string& abortMessage() const { return abort_message_; }

    long long pointsSolved() const { return points_solved_; }
    long long newtonIterations() const { return newton_iterations_; }
    long long predictorHits() const { return predictor_hits_; }
    int retrySolves() const { return retry_solves_; }
    int failedPoints() const { return failed_points_; }
    int maxPointIterations() const { return max_point_iterations_; }

private:
    /**
     * @brief Solves the point at inner-level value v. On success it
     *        publishes the solution to the circuit and advances the chain.
     */
    NLSolveResult solveAt(double v, const std::vector<double>& x0) {
        // primary_ is the stepped variable: the inner one for a nested
        // sweep, the swept one itself for a flat sweep (the driver picks it
        // that way). secondary_, when present, is held by solveUnit().
        primary_.set(v);
        const NLSolveResult nl = system_.solveNonlinear(/*time=*/0.0,
                                                        /*q_coeff=*/0.0,
                                                        /*use_history=*/false,
                                                        x0);
        if (nl.converged) {
            system_.extractSolution();
            x_ante_ = x_last_;
            ante_v_ = last_v_;
            x_last_ = system_.solution();
            last_v_ = v;
            chain_ = true;
            last_iters_ = nl.total_iterations;
            newton_iterations_ += nl.total_iterations;
            max_point_iterations_ =
                std::max(max_point_iterations_, nl.total_iterations);
        }
        return nl;
    }

    DaeSystem& system_;
    Circuit& circuit_;
    SweepTarget& primary_;
    SweepTarget* secondary_;
    EngineConfig cfg_;

    // Continuation state (last two accepted points of the chain).
    std::vector<double> x_last_, x_ante_, seed_, cold_, row_;
    double last_v_ = 0.0, ante_v_ = 0.0;
    bool   chain_ = false;
    int    last_iters_ = 0;

    long long points_solved_ = 0;
    long long newton_iterations_ = 0;
    long long predictor_hits_ = 0;
    int retry_solves_ = 0;
    int failed_points_ = 0;
    int max_point_iterations_ = 0;

    bool aborted_ = false;
    double abort_value_ = 0.0, abort_outer_ = 0.0;
    size_t abort_index_ = 0;
    std::string abort_message_;
};

} // namespace

DCSweep::DCSweep(Circuit& circuit)
    : circuit_(circuit)
{
}

DCSweep::~DCSweep() = default;

std::unique_ptr<SweepTarget> DCSweep::makeTarget(Circuit& circuit,
                                                 const std::string& name,
                                                 std::string& error) {
    Device* raw = circuit.getDevice(name);
    if (raw != nullptr) {
        if (auto* vs = dynamic_cast<VoltageSource*>(raw)) {
            return std::unique_ptr<SweepTarget>(new VSourceTarget(vs));
        }
        if (auto* is = dynamic_cast<CurrentSource*>(raw)) {
            return std::unique_ptr<SweepTarget>(new ISourceTarget(is));
        }
        if (auto* r = dynamic_cast<Resistor*>(raw)) {
            return std::unique_ptr<SweepTarget>(new ResistorTarget(r));
        }
        error = "'" + name + "' is not sweepable (expected an independent "
                "voltage source, current source or resistor)";
        return nullptr;
    }
    // Xyce: the swept variable may also be TEMP (circuit temperature, in
    // Celsius). A device actually named TEMP wins — a real name is never
    // shadowed by the keyword.
    std::string upper = name;
    for (char& c : upper) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    if (upper == "TEMP") {
        return std::unique_ptr<SweepTarget>(new TemperatureTarget(&circuit));
    }
    error = "no device named '" + name + "' in the circuit (and it is not TEMP)";
    return nullptr;
}

void DCSweep::setCsvStream(std::ostream* os, std::vector<std::string> labels) {
    csv_ = os;
    csv_labels_ = std::move(labels);
    csv_header_written_ = false;
}

std::vector<double> DCSweep::buildGrid(const DCSweepSpec& spec) const {
    std::vector<double> grid;
    if (spec.variation == DCSweepVariation::LIN) {
        const double step = spec.step;
        if (step == 0.0) {                 // degenerate: a single point
            grid.push_back(spec.start);
            return grid;
        }
        grid.reserve(64);
        const double eps = std::abs(spec.stop) * 1e-12 + 1e-300;
        double v = spec.start;
        if (step > 0.0) {
            while (v <= spec.stop + eps) {
                grid.push_back(v);
                if (static_cast<int>(grid.size()) >= kMaxPoints_) break;
                const double next = v + step;
                if (next <= v) break;       // step underflow guard
                v = next;
            }
        } else {
            while (v >= spec.stop - eps) {
                grid.push_back(v);
                if (static_cast<int>(grid.size()) >= kMaxPoints_) break;
                const double next = v + step;
                if (next >= v) break;
                v = next;
            }
        }
        return grid;
    }

    // DEC/OCT logarithmic grids. Xyce requires start > 0 and stop > start;
    // otherwise it runs a single point at start (warning issued by run()).
    const int points = std::max(1, static_cast<int>(std::llround(spec.step)));
    if (spec.start <= 0.0 || spec.stop <= spec.start) {
        grid.push_back(spec.start);
        return grid;
    }
    const double base = (spec.variation == DCSweepVariation::DEC) ? 10.0 : 2.0;
    const double decades = std::log(spec.stop / spec.start) / std::log(base);
    const long long total =
        static_cast<long long>(std::ceil(decades * points)) + 1;
    grid.reserve(static_cast<size_t>(std::min<long long>(total, kMaxPoints_)));
    for (long long k = 0; k < total; ++k) {
        if (static_cast<int>(grid.size()) >= kMaxPoints_) break;
        const double v = spec.start * std::pow(base, static_cast<double>(k) / points);
        if (v > spec.stop * (1.0 + 1e-9)) break;
        grid.push_back(v);
    }
    // Pin the exact stop (avoids a short last interval from fp rounding).
    if (!grid.empty() && grid.back() < spec.stop * (1.0 - 1e-9)) {
        grid.push_back(spec.stop);
    }
    return grid;
}

DCSweepResult DCSweep::run(bool use_mpi) {
    DCSweepResult result;
    result.sweep_source = spec_.source;
    result.variation    = variationName(spec_.variation);
    result.nested       = spec2_.valid();
    if (result.nested) {
        result.sweep_source2 = spec2_.source;
        result.variation2    = variationName(spec2_.variation);
    }

    auto fail = [&](const std::string& msg) {
        result.success = false;
        result.message = msg;
        return result;
    };

    if (!spec_.valid()) return fail("no DC sweep variable set");

    std::string error_msg;
    if (!circuit_.validate(error_msg)) {
        return fail("Circuit validation failed: " + error_msg);
    }

    std::unique_ptr<SweepTarget> outer = makeTarget(circuit_, spec_.source, error_msg);
    if (!outer) {
        return fail("DC sweep variable error: " + error_msg);
    }

    std::unique_ptr<SweepTarget> inner;
    if (result.nested) {
        if (spec2_.source == spec_.source) {
            return fail("nested .DC sweep repeats the same variable '"
                        + spec_.source + "'");
        }
        inner = makeTarget(circuit_, spec2_.source, error_msg);
        if (!inner) {
            return fail("nested DC sweep variable error: " + error_msg);
        }
    }

    const std::vector<double> grid_outer = buildGrid(spec_);
    if (grid_outer.empty()) return fail("DC sweep grid is empty");
    std::vector<double> grid_inner;
    if (result.nested) {
        grid_inner = buildGrid(spec2_);
        if (grid_inner.empty()) return fail("nested DC sweep grid is empty");
    }
    const size_t n_inner = result.nested ? grid_inner.size() : 1;
    const size_t total_points = grid_outer.size() * n_inner;

    // Chain units: one per point (flat) or one per outer row (nested).
    const size_t n_units = grid_outer.size();
    const size_t units_per = n_inner;

    if (verbose_) {
        std::cout << "[DCSweep] " << spec_.source << " (" << outer->kind()
                  << ") " << grid_outer.size() << " points";
        if (result.nested) {
            std::cout << " x " << spec2_.source << " (" << inner->kind() << ") "
                      << grid_inner.size();
        }
        std::cout << " = " << total_points << " operating points"
                  << (continue_on_fail_ ? " [NOFAIL]" : "") << std::endl;
        if (result.nested && grid_inner.size() == 1
            && spec2_.variation != DCSweepVariation::LIN) {
            std::cout << "[DCSweep] warning: nested range degenerate, "
                         "single point at start" << std::endl;
        }
    }

    // Netlist values are restored on every exit path (error included).
    const double saved_outer = outer->get();
    const double saved_inner = inner ? inner->get() : 0.0;
    ScopeGuard restore;
    restore.fn = [&] {
        outer->set(saved_outer);
        if (inner) inner->set(saved_inner);
    };

    // ---- ONE setup for the whole sweep -------------------------------------
    // Across points only source values move: pattern, CSR slots, StateStore
    // layout and the MPI distribution plan are invariant, so they are built
    // once instead of once per point.
    DaeSystem system(circuit_);
    // Point distribution replicates the circuit per rank, so the inner
    // device/row distribution must be off (each rank needs the whole
    // solution to report its own points).
    const bool distribute_points =
        point_distribution_ && use_mpi
        && MPIManager::getInstance().getSize() > 1
        && total_points >= static_cast<size_t>(2)
                            * static_cast<size_t>(MPIManager::getInstance().getSize());
    system.setMPIManager((use_mpi && !distribute_points)
                         ? &MPIManager::getInstance() : nullptr);
    system.setup();
    system.setMaxIterations(max_iterations_);
    system.setRelTol(1e-3);
    system.setAbsTol(std::max(1e-12, tolerance_ * 1e-3));
    system.setGminSchedule(1e-3, 1e-12);
    system.setSolverType(solver_type_);
    system.setVerbose(verbose_);
    system.beginSolve(/*initial_solve=*/true);
    system.resetSolveStats();

    const auto t0 = std::chrono::high_resolution_clock::now();

    EngineConfig cfg;
    cfg.predictor = predictor_;
    cfg.continue_on_fail = continue_on_fail_;
    cfg.verbose = verbose_;
    cfg.max_retries = max_retries_;

    SweepTarget* primary = inner.get() ? inner.get() : outer.get();
    SweepTarget* secondary = inner.get() ? outer.get() : nullptr;

    // In-memory retention stride (RAM is capped; the CSV keeps every point).
    size_t stride = 1;
    if (max_retained_ > 0 && total_points > max_retained_) {
        stride = (total_points + max_retained_ - 1) / max_retained_;
    }
    const auto retain = [this, stride, total_points](size_t index) -> bool {
        return max_retained_ == 0 || (index % stride == 0)
               || (index + 1 == total_points);
    };

    // Full precision for the streamed rows (default ostream precision would
    // silently truncate the sweep data to 6 significant digits).
    if (csv_) {
        csv_saved_precision_ = csv_->precision(17);
    }

    const std::string inner_name = inner ? spec2_.source : spec_.source;

    auto writeHeader = [&] {
        if (!csv_ || csv_header_written_) return;
        *csv_ << spec_.source;
        if (result.nested) *csv_ << "," << spec2_.source;
        for (const auto& label : csv_labels_) *csv_ << "," << label;
        *csv_ << "\n";
        csv_header_written_ = true;
    };

    // Reports one solved point: CSV row + (decimated) in-memory sample.
    auto report = [&](const PointRow& r) {
        if (csv_) {
            writeHeader();
            if (result.nested) *csv_ << r.outer << "," << r.value;
            else               *csv_ << r.value;
            for (double x : r.nodes) *csv_ << "," << x;
            *csv_ << "\n";
        }
        if (!retain(r.index)) return;
        DCSweepPoint point;
        point.sweep_value = r.value;
        point.outer_value = r.outer;
        point.node_voltages = r.nodes;
        point.branch_currents = r.currents;
        point.newton_iterations = r.iters;
        point.converged = r.converged;
        result.points.push_back(std::move(point));
    };

    // Outcome of the distributed path, applied in the common epilogue below.
    bool dist_success = true;
    double dist_fail_value = 0.0, dist_fail_outer = 0.0;
    std::vector<double> dist_decision;

    if (!distribute_points) {
        // ================= sequential chain (default) ======================
        // Every rank walks the same units in the same order with the same
        // seeds: identical numbers everywhere, no extra communication.
        SweepEngine engine(system, circuit_, *primary, secondary, cfg);
        bool aborted = false;
        for (size_t u = 0; u < n_units && !aborted; ++u) {
            std::vector<PointRow> rows;
            rows.reserve(units_per);
            aborted = !engine.solveUnit(grid_outer, grid_inner, u,
                                       u * units_per, rows);
            for (auto& r : rows) {
                if (retain(r.index)) engine.addCurrents(r);
                report(r);
            }
            if (verbose_ && (u % 100 == 0 || u + 1 == n_units)) {
                std::cout << "[DCSweep] point " << ((u + 1) * units_per)
                          << "/" << total_points << " " << inner_name << "="
                          << grid_outer[u] << " (chain unit " << (u + 1)
                          << "/" << n_units << ")" << std::endl;
            }
            if (progress_ && (u % 100 == 0 || u + 1 == n_units)) {
                progress_((u + 1) * units_per, total_points, grid_outer[u]);
            }
        }

        if (aborted && !continue_on_fail_) {
            return fail("DC sweep failed at " + inner_name + "="
                        + std::to_string(engine.abortValue())
                        + (result.nested
                           ? " (" + spec_.source + "="
                             + std::to_string(engine.abortOuter()) + ")" : "")
                        + " (point " + std::to_string(engine.abortIndex() + 1)
                        + "/" + std::to_string(total_points) + "): "
                        + engine.abortMessage());
        }
        if (verbose_ && continue_on_fail_ && engine.failedPoints() > 0) {
            std::cout << "[DCSweep] warning: " << engine.failedPoints()
                      << " point(s) without solution (NOFAIL)" << std::endl;
        }

        result.stats.points_solved = engine.pointsSolved();
        result.stats.newton_iterations = engine.newtonIterations();
        result.stats.predictor_hits = engine.predictorHits();
        result.stats.retry_solves = engine.retrySolves();
        result.stats.failed_points = engine.failedPoints();
        result.stats.max_point_iterations = engine.maxPointIterations();
        result.stats.linear_solves = system.solveStats().linearSolves;
        result.stats.ranks_used = 1;
    } else {
        // ================= points distributed across ranks =================
        // Each rank owns a contiguous block of chain units (contiguity keeps
        // the local continuation chain dense) and runs its own replica of
        // the circuit. Points are then exchanged one by one so the master
        // writes the CSV in sweep order.
        MPIManager& mpi = MPIManager::getInstance();
        const int rank = mpi.getRank();
        const int size = mpi.getSize();
        int unit_begin = 0, unit_end = 0;
        mpi.blockRangeForRank(static_cast<int>(n_units), rank,
                              unit_begin, unit_end);

        SweepEngine engine(system, circuit_, *primary, secondary, cfg);
        std::vector<PointRow> mine;                 // this rank's solved points
        for (int u = unit_begin; u < unit_end; ++u) {
            engine.solveUnit(grid_outer, grid_inner, static_cast<size_t>(u),
                             static_cast<size_t>(u) * units_per, mine);
            if (engine.aborted() && !continue_on_fail_) break;  // stay collective
        }
        for (PointRow& r : mine) {
            if (retain(r.index)) engine.addCurrents(r);
        }

        // ---- exchange in BATCHES -------------------------------------------
        // Ranks own CONTIGUOUS unit blocks, so a global point index maps to a
        // local offset arithmetically (O(1), no search). Points travel in
        // chunks, so the sweep costs a handful of collectives instead of one
        // per point — that batching is what makes the mode competitive.
        //
        // Payload per chunk: [n_points, then per point:
        //                    index, value, outer, iters, converged,
        //                    <node voltages>, n_currents, <currents>]
        const int nnodes = std::max(0, system.numNodeRows());
        const size_t first_point = static_cast<size_t>(unit_begin) * units_per;
        const size_t last_point  = static_cast<size_t>(unit_end) * units_per;
        constexpr size_t kExchangeChunk = 256;

        std::vector<double> local;
        std::vector<double> gathered;
        size_t emitted = 0;
        for (size_t chunk = 0; chunk < total_points; chunk += kExchangeChunk) {
            const size_t hi = std::min(total_points, chunk + kExchangeChunk);
            local.clear();
            const size_t lo = std::max(chunk, first_point);
            const size_t bh = std::min(hi, last_point);
            if (bh > lo) {
                local.push_back(static_cast<double>(bh - lo));
                for (size_t i = lo; i < bh; ++i) {
                    const PointRow& r = mine[i - first_point];
                    local.push_back(static_cast<double>(r.index));
                    local.push_back(r.value);
                    local.push_back(r.outer);
                    local.push_back(static_cast<double>(r.iters));
                    local.push_back(r.converged ? 1.0 : 0.0);
                    local.insert(local.end(), r.nodes.begin(), r.nodes.end());
                    local.push_back(static_cast<double>(r.currents.size()));
                    local.insert(local.end(), r.currents.begin(),
                                 r.currents.end());
                }
            }
            mpi.gather(local, gathered);

            if (!mpi.isMaster()) continue;
            // Arrival order == ascending global index (rank order == unit
            // order for contiguous blocks), so the CSV keeps sweep order.
            size_t off = 0;
            while (off < gathered.size()) {
                const size_t n = static_cast<size_t>(gathered[off++]);
                for (size_t k = 0; k < n && off + 5.0 <= gathered.size(); ++k) {
                    PointRow r;
                    r.index     = static_cast<size_t>(gathered[off + 0]);
                    r.value     = gathered[off + 1];
                    r.outer     = gathered[off + 2];
                    r.iters     = static_cast<int>(gathered[off + 3]);
                    r.converged = gathered[off + 4] != 0.0;
                    off += 5;
                    const size_t avail = gathered.size() - off;
                    const size_t take = std::min<size_t>(
                        static_cast<size_t>(nnodes), avail);
                    r.nodes.assign(gathered.begin() + static_cast<long>(off),
                                   gathered.begin() + static_cast<long>(off + take));
                    off += take;
                    if (avail > take && off < gathered.size()) {
                        const size_t nc = static_cast<size_t>(gathered[off++]);
                        const size_t take2 =
                            std::min<size_t>(nc, gathered.size() - off);
                        r.currents.assign(
                            gathered.begin() + static_cast<long>(off),
                            gathered.begin() + static_cast<long>(off + take2));
                        off += take2;
                    }
                    if (r.index != emitted) continue;   // defensive: no dupes
                    report(r);
                    ++emitted;
                    if (verbose_ && (emitted % 100 == 0 || emitted == total_points)) {
                        std::cout << "[DCSweep] point " << emitted << "/"
                                  << total_points << " " << inner_name << "="
                                  << r.value << std::endl;
                    }
                    if (progress_ && (emitted % 100 == 0
                                      || emitted == total_points)) {
                        progress_(emitted, total_points, r.value);
                    }
                }
            }
        }

        // ---- reduce statistics and failure information -------------------
        const double local_stats[7] = {
            static_cast<double>(engine.pointsSolved()),
            static_cast<double>(engine.newtonIterations()),
            static_cast<double>(engine.predictorHits()),
            static_cast<double>(engine.retrySolves()),
            static_cast<double>(engine.failedPoints()),
            static_cast<double>(engine.maxPointIterations()),
            static_cast<double>(system.solveStats().linearSolves)
        };
        std::vector<double> all_stats;
        mpi.gather(std::vector<double>(local_stats, local_stats + 7), all_stats);
        const double fail_info[4] = {
            (engine.aborted() || engine.failedPoints() > 0) ? 1.0 : 0.0,
            engine.abortValue(),
            engine.abortOuter(),
            static_cast<double>(engine.abortIndex())
        };
        std::vector<double> all_fail;
        mpi.gather(std::vector<double>(fail_info, fail_info + 4), all_fail);

        // Decision vector: the master owns the outcome and then tells every
        // rank, so ALL ranks report the same success/message (main returns a
        // PER-RANK exit status). Points and the CSV stay master-only.
        std::vector<double> decision(12, 0.0);
        if (mpi.isMaster()) {
            long long sum[5] = {0, 0, 0, 0, 0};
            int max_iters = 0;
            long long linear = 0;
            for (size_t k = 0; k + 6 < all_stats.size(); k += 7) {
                sum[0] += static_cast<long long>(all_stats[k + 0]);
                sum[1] += static_cast<long long>(all_stats[k + 1]);
                sum[2] += static_cast<long long>(all_stats[k + 2]);
                sum[3] += static_cast<long long>(all_stats[k + 3]);
                sum[4] += static_cast<long long>(all_stats[k + 4]);
                max_iters = std::max(max_iters,
                                     static_cast<int>(all_stats[k + 5]));
                linear  += static_cast<long long>(all_stats[k + 6]);
            }
            bool any_fail = false;
            for (size_t k = 0; k + 3 < all_fail.size(); k += 4) {
                if (all_fail[k] == 0.0) continue;
                if (!any_fail) {     // lowest rank owns the earliest point
                    decision[9]  = all_fail[k + 1];
                    decision[10] = all_fail[k + 2];
                    decision[11] = all_fail[k + 3] + 1.0;
                }
                any_fail = true;
            }
            decision[0] = (any_fail && !continue_on_fail_) ? 0.0 : 1.0;
            decision[1] = static_cast<double>(sum[0]);
            decision[2] = static_cast<double>(sum[1]);
            decision[3] = static_cast<double>(sum[2]);
            decision[4] = static_cast<double>(sum[3]);
            decision[5] = static_cast<double>(sum[4]);
            decision[6] = static_cast<double>(max_iters);
            decision[7] = static_cast<double>(linear);
            decision[8] = static_cast<double>(size);
        }
        mpi.broadcast(decision);

        dist_success = decision[0] != 0.0;
        dist_decision = decision;
        dist_fail_value = decision[9];
        dist_fail_outer = decision[10];
    }

    // ---- common epilogue ---------------------------------------------------
    if (distribute_points) {
        result.stats.points_solved =
            static_cast<long long>(dist_decision[1]);
        result.stats.newton_iterations =
            static_cast<long long>(dist_decision[2]);
        result.stats.predictor_hits =
            static_cast<long long>(dist_decision[3]);
        result.stats.retry_solves =
            static_cast<int>(dist_decision[4]);
        result.stats.failed_points =
            static_cast<int>(dist_decision[5]);
        result.stats.max_point_iterations =
            static_cast<int>(dist_decision[6]);
        result.stats.linear_solves =
            static_cast<long long>(dist_decision[7]);
        result.stats.ranks_used = static_cast<int>(dist_decision[8]);
        if (!dist_success) {
            return fail("DC sweep failed at " + inner_name + "="
                        + std::to_string(dist_fail_value)
                        + (result.nested
                           ? " (" + spec_.source + "="
                             + std::to_string(dist_fail_outer) + ")" : "")
                        + " (point "
                        + std::to_string(static_cast<long long>(dist_decision[11]))
                        + "/" + std::to_string(total_points)
                        + ", distributed across ranks)");
        }
    }

    result.stats.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();

    if (csv_) {
        csv_->precision(csv_saved_precision_);
        csv_->flush();
    }

    result.success = true;
    result.message = "DC sweep converged (" +
        std::to_string(result.stats.points_solved) + " points";
    if (result.stats.failed_points > 0) {
        result.message += ", " + std::to_string(result.stats.failed_points) +
                          " without solution";
    }
    result.message += ")";
    if (result.stats.retry_solves > 0) {
        result.message += " [" + std::to_string(result.stats.retry_solves) +
                          " source-step retries]";
    }
    if (result.stats.ranks_used > 1) {
        result.message += " [points distributed over " +
                          std::to_string(result.stats.ranks_used) + " ranks]";
    }
    if (verbose_) {
        std::cout << "[DCSweep] " << result.message << ", "
                  << result.stats.newton_iterations << " Newton iterations"
                  << " (" << result.stats.linear_solves << " linear solves), "
                  << "predictor used on " << result.stats.predictor_hits
                  << " points, " << std::fixed << std::setprecision(3)
                  << result.stats.elapsed_seconds << " s" << std::endl;
    }
    return result;
}

} // namespace e_XSpice