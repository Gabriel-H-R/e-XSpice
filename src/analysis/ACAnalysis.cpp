/**
 * @author Gabriel Henrique Silva
 * @file ACAnalysis.cpp
 * @brief Small-signal AC analysis: operating point once, complex solve per
 *        frequency.
 */

#include "ACAnalysis.h"

#include "../core/Device.h"
#include "../parallel/MPIManager.h"
#include "../solver/ComplexSparseLU.h"
#include "SmallSignalLinearization.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace e_XSpice {

namespace {
constexpr double kPi2 = 6.28318530717958647692;   // 2*pi
constexpr double kDegPerRad = 57.29577951308232088;
}  // namespace

ACAnalysis::ACAnalysis(Circuit& circuit) : circuit_(circuit) {}

ACAnalysis::~ACAnalysis() {
    if (csv_) csv_->precision(csv_saved_precision_);
}

// ===============================================================
// Frequency grid
// ===============================================================

std::vector<double> ACAnalysis::buildGrid(ACVariation variation, int points,
                                          double fstart, double fstop) {
    std::vector<double> grid;
    if (points < 1 || !(fstart > 0.0) || !(fstop > 0.0)) return grid;

    if (variation == ACVariation::LIN) {
        const int n = std::max(1, points);
        grid.reserve(static_cast<size_t>(n));
        if (n == 1) {
            grid.push_back(fstart);
        } else {
            const double step = (fstop - fstart) / static_cast<double>(n - 1);
            for (int i = 0; i < n; ++i) {
                grid.push_back(fstart + step * i);
            }
            grid.back() = fstop;     // exact endpoint, no accumulation drift
        }
        return grid;
    }

    // DEC / OCT: log-spaced, both endpoints included, `points` per decade or
    // octave (same convention as the .DC sweep and Xyce).
    const double mult = (variation == ACVariation::DEC) ? 10.0 : 2.0;
    const double lo = std::log(fstart);
    const double hi = std::log(fstop);
    const double span = (hi - lo) / std::log(mult);
    if (!(span > 0.0)) {
        grid.push_back(fstart);
        return grid;
    }
    const int n = std::max(1, points);
    const long long steps = static_cast<long long>(std::floor(span * n + 0.5));
    grid.reserve(static_cast<size_t>(steps) + 1);
    for (long long i = 0; i <= steps; ++i) {
        const double f = std::exp(lo + static_cast<double>(i) / n
                                  * std::log(mult));
        if (f > fstop) break;
        grid.push_back(f);
    }
    // The last computed point can land a rounding step short of fstop; snap it
    // instead of appending a duplicate (and pin the endpoint exactly, so the
    // CSV is reproducible run to run).
    if (!grid.empty()) {
        const double rel = std::fabs(grid.back() - fstop) / fstop;
        if (rel <= 1e-12) {
            grid.back() = fstop;
        } else if (grid.back() < fstop) {
            grid.push_back(fstop);
        }
    } else {
        grid.push_back(fstart);
    }
    return grid;
}

void ACAnalysis::setCsvStream(std::ostream* os, std::vector<std::string> labels) {
    csv_ = os;
    csv_labels_ = std::move(labels);
    csv_header_written_ = false;
}

// ===============================================================
// Run
// ===============================================================

ACResult ACAnalysis::run(bool use_mpi) {
    ACResult result;
    result.variation = (variation_ == ACVariation::LIN)  ? "LIN"
                     : (variation_ == ACVariation::DEC)  ? "DEC" : "OCT";

    const std::vector<double> grid =
        buildGrid(variation_, points_, fstart_, fstop_);
    if (grid.empty()) {
        result.message = "invalid frequency sweep specification";
        return result;
    }
    const size_t nfreq = grid.size();
    if (nfreq > static_cast<size_t>(kMaxFrequencies_)) {
        result.message = "frequency grid too large (max " +
                         std::to_string(kMaxFrequencies_) + " points)";
        return result;
    }

    MPIManager& mpi = MPIManager::getInstance();
    const int rank = mpi.getRank();
    const int world = mpi.getSize();

    // Frequencies are distributed only when it can actually pay off: with
    // fewer frequencies than ranks the exchange would dominate.
    const bool distribute = point_distribution_ && use_mpi && world > 1 &&
                            nfreq >= static_cast<size_t>(2) *
                                         static_cast<size_t>(world);

    if (verbose_ && rank == 0) {
        std::cout << "[ACAnalysis] " << result.variation << " sweep "
                  << nfreq << " frequencies, " << fstart_ << " Hz .. "
                  << fstop_ << " Hz" << std::endl;
    }

    // ---- operating point + linearization --------------------------------
    // The AC analysis owns its operating point (SPICE/Xyce always compute one
    // before the sweep). SmallSignalLinearization does the whole preparation —
    // OP solve, (G, C) snapshot, ground elimination, pivot-safe ordering — and
    // is the SAME object the noise analysis uses, so both analyses provably
    // see the same frequency system.
    const auto t0 = std::chrono::high_resolution_clock::now();

    SmallSignalLinearization lin;
    {
        SmallSignalLinearization::Options opt;
        opt.max_iterations = max_iterations_;
        opt.abs_tol = std::max(1e-12, tolerance_ * 1e-3);
        opt.verbose = verbose_ && rank == 0;
        opt.use_distributed_layout = (use_mpi && !distribute);
        if (!lin.prepare(circuit_, opt)) {
            result.success = false;
            result.message = lin.message();
            result.stats.elapsed_seconds =
                std::chrono::duration<double>(
                    std::chrono::high_resolution_clock::now() - t0).count();
            if (rank == 0) std::cerr << result.message << std::endl;
            return result;
        }
    }
    result.stats.newton_iterations = lin.newtonIterations();

    if (verbose_ && rank == 0) {
        std::cout << "[ACAnalysis] operating point in "
                  << result.stats.newton_iterations << " Newton iterations"
                  << std::endl;
    }

    const int n_full = lin.fullSize();
    const int num_node_rows = lin.numNodeRows();
    const int nn = lin.size();

    std::vector<int> ac_of_row(static_cast<size_t>(n_full), -1);
    std::vector<int> row_of_ac(static_cast<size_t>(n_full), -1);
    for (int r = 0; r < n_full; ++r) ac_of_row[static_cast<size_t>(r)] = lin.acRow(r);
    for (int k = 0; k < nn; ++k) row_of_ac[static_cast<size_t>(k)] = lin.fullRow(k);

    const std::vector<int>&    ac_row_ptr = lin.rowPointers();
    const std::vector<int>&    ac_col_idx = lin.colIndices();
    const std::vector<double>& g_vals = lin.gValues();
    const std::vector<double>& c_vals = lin.cValues();
    const StaticRowPermutation& perm = lin.permutation();

    // Branch-current labels: extra unknowns are assigned in device order (row
    // num_node_rows + k), so I(V1) / I(L1) ... line up with this list.
    std::vector<std::string> branch_labels;
    for (const auto& d : circuit_.getDevices()) {
        const int extra = std::max(0, d->numExtraUnknowns());
        for (int k = 0; k < extra; ++k) {
            branch_labels.push_back(d->getName());
        }
    }
    const int num_branch = std::max(0, n_full - num_node_rows);

    // Node phasors are reported per NON-GROUND node, in Circuit::getNodes()
    // order — the same order the CSV labels and the log use.
    std::vector<int> node_ac_row;
    node_ac_row.reserve(static_cast<size_t>(num_node_rows));
    for (const auto& n : circuit_.getNodes()) {
        if (n->isGround()) continue;
        const int id = n->getId();
        const auto& node_map = lin.system().nodeMap();
        const int row = (id >= 0 && id < static_cast<int>(node_map.size()))
                            ? node_map[static_cast<size_t>(id)] : -1;
        node_ac_row.push_back(ac_of_row[static_cast<size_t>(row < 0 ? 0 : row)]);
    }
    const int n_nodes_reported = static_cast<int>(node_ac_row.size());

    // AC row of every extra unknown (branch current).
    std::vector<int> branch_ac_row(static_cast<size_t>(num_branch), -1);
    for (int k = 0; k < num_branch; ++k) {
        branch_ac_row[static_cast<size_t>(k)] =
            ac_of_row[static_cast<size_t>(num_node_rows + k)];
    }

    // ---- excitation -----------------------------------------------------
    // The excitation is independent of frequency, so it is assembled ONCE (a
    // source phasor is frequency independent by definition). Every rank folds
    // the full device list into its own replicated buffers: no communication,
    // and identical on every rank by construction. The rows/columns of the
    // eliminated ground are dropped here too.
    std::vector<double> b_re(static_cast<size_t>(nn), 0.0);
    std::vector<double> b_im(static_cast<size_t>(nn), 0.0);
    {
        std::vector<double> full_re(static_cast<size_t>(n_full), 0.0);
        std::vector<double> full_im(static_cast<size_t>(n_full), 0.0);
        for (const auto& d : circuit_.getDevices()) {
            d->loadACExcitation(1.0, full_re, full_im);
        }
        b_re = lin.permuteRhs(full_re);
        b_im = lin.permuteRhs(full_im);
    }

    // ---- retention ------------------------------------------------------
    size_t stride = 1;
    if (max_retained_ > 0 && nfreq > max_retained_) {
        stride = (nfreq + max_retained_ - 1) / max_retained_;
    }
    const auto retain = [this, stride, nfreq](size_t index) -> bool {
        return max_retained_ == 0 || (index % stride == 0)
               || (index + 1 == nfreq);
    };

    if (csv_) csv_saved_precision_ = csv_->precision(17);

    if (std::getenv("EXSPICE_DEBUG_AC")) {
        std::fprintf(stderr, "[AC] n_full=%d gnd_row=%d nn=%d perm=", n_full,
                     lin.groundRow(), nn);
        for (int k = 0; k < nn; ++k) {
            std::fprintf(stderr, " %d", perm.row_of_pos[static_cast<size_t>(k)]);
        }
        std::fprintf(stderr, "\n[AC] node rows:");
        for (int a : node_ac_row) std::fprintf(stderr, " %d", a);
        std::fprintf(stderr, "\n[AC] branch rows:");
        for (int a : branch_ac_row) std::fprintf(stderr, " %d", a);
        std::fprintf(stderr, "\n");
    }

    auto writeHeader = [&] {
        if (!csv_ || csv_header_written_) return;
        *csv_ << "frequency";
        for (size_t i = 0; i < csv_labels_.size(); ++i) {
            *csv_ << "," << csv_labels_[i] << "_mag"
                  << "," << csv_labels_[i] << "_phase";
        }
        for (const auto& l : branch_labels) {
            *csv_ << "," << l << "_imag_mag" << "," << l << "_imag_phase";
        }
        *csv_ << "\n";
        csv_header_written_ = true;
    };

    // One solved frequency, held as raw phasors (the magnitude/phase post
    // processing happens at reporting time).
    struct Row {
        size_t index = 0;
        bool  valid = false;
        std::vector<double> re;   // length nn
        std::vector<double> im;
    };

    int lo_i = 0, hi_i = static_cast<int>(nfreq);
    if (distribute) mpi.blockRangeForRank(static_cast<int>(nfreq), rank, lo_i, hi_i);

    long long flops_total = 0;
    int failures_total = 0;

    // ---- per-thread factor slots ---------------------------------------
    // Frequencies are independent, so they are solved CONCURRENTLY on the
    // replicated (G, C). Each thread owns one ComplexSparseLU (its own factor
    // arrays) plus its own solution buffers: no sharing, no lock, and the
    // result does not depend on how the frequencies were grouped.
    int nthreads = 1;
    #ifdef _OPENMP
    nthreads = std::max(1, omp_get_max_threads());
    #endif
    const int n_local = hi_i - lo_i;
    const bool parallel = (nthreads > 1) && (n_local > 1);
    const int nslots = parallel ? std::min(nthreads, n_local) : 1;

    // Per-node peak accumulators (filled by report()).
    result.stats.node_peak_mag.assign(static_cast<size_t>(n_nodes_reported), 0.0);
    result.stats.node_peak_db.assign(static_cast<size_t>(n_nodes_reported), 0.0);
    result.stats.node_peak_freq.assign(static_cast<size_t>(n_nodes_reported), 0.0);

    std::vector<ComplexSparseLU> lus(static_cast<size_t>(nslots));
    for (int s = 0; s < nslots; ++s) {
        // The symbolic factorization runs once per slot; every frequency then
        // only refills values on this static pattern (the PERMUTED one — it
        // must be the same pattern the value arrays follow).
        lus[static_cast<size_t>(s)].setPattern(ac_row_ptr, ac_col_idx);
    }

    // Frequencies are processed in BLOCKS so RAM stays bounded: a block holds
    // at most ~2^20 phasor doubles regardless of grid length and system size.
    const size_t payload_cap = static_cast<size_t>(1) << 20;
    const size_t block = std::max<size_t>(
        1, std::min<size_t>(4096, payload_cap / std::max<size_t>(2, static_cast<size_t>(nn))));

    // ---- report ---------------------------------------------------------
    // The reporting path is identical for the local and the distributed mode
    // (the master receives the rows and fills the same ACDataPoint), so the
    // CSV and the retained sample cannot diverge between the two.
    auto report = [&](const Row& r) {
        const double f = grid[r.index];
        ACDataPoint p;
        p.frequency = f;
        p.valid = r.valid;
        p.mag.resize(static_cast<size_t>(n_nodes_reported));
        p.phase_deg.resize(static_cast<size_t>(n_nodes_reported));
        p.db.resize(static_cast<size_t>(n_nodes_reported));
        p.branch_mag.resize(static_cast<size_t>(num_branch));
        p.branch_phase_deg.resize(static_cast<size_t>(num_branch));
        for (int i = 0; i < n_nodes_reported; ++i) {
            const int a = node_ac_row[static_cast<size_t>(i)];
            const double re = (a >= 0) ? r.re[static_cast<size_t>(a)] : 0.0;
            const double im = (a >= 0) ? r.im[static_cast<size_t>(a)] : 0.0;
            const double m = std::hypot(re, im);
            p.mag[static_cast<size_t>(i)] = m;
            p.phase_deg[static_cast<size_t>(i)] = std::atan2(im, re) * kDegPerRad;
            p.db[static_cast<size_t>(i)] = 20.0 * std::log10(std::max(m, 1e-300));
        }
        for (int k = 0; k < num_branch; ++k) {
            const int a = branch_ac_row[static_cast<size_t>(k)];
            const double re = (a >= 0) ? r.re[static_cast<size_t>(a)] : 0.0;
            const double im = (a >= 0) ? r.im[static_cast<size_t>(a)] : 0.0;
            const double m = std::hypot(re, im);
            p.branch_mag[static_cast<size_t>(k)] = m;
            p.branch_phase_deg[static_cast<size_t>(k)] =
                std::atan2(im, re) * kDegPerRad;
        }
        writeHeader();
        // Resonance summary: tracked over EVERY solved frequency (not the
        // decimated sample), so a long sweep still reports its true peak.
        for (int i = 0; i < n_nodes_reported; ++i) {
            const size_t ui = static_cast<size_t>(i);
            const double m = p.mag[ui];
            if (m > result.stats.node_peak_mag[ui]) {
                result.stats.node_peak_mag[ui] = m;
                result.stats.node_peak_db[ui] = p.db[ui];
                result.stats.node_peak_freq[ui] = f;
            }
            if (m > result.stats.peak_mag) {
                result.stats.peak_mag = m;
                result.stats.peak_db = p.db[ui];
                result.stats.peak_frequency = f;
                result.stats.peak_node = i;
            }
        }
        if (csv_) {
            *csv_ << f;
            for (int i = 0; i < n_nodes_reported; ++i) {
                *csv_ << "," << p.mag[static_cast<size_t>(i)]
                      << "," << p.phase_deg[static_cast<size_t>(i)];
            }
            for (int k = 0; k < num_branch; ++k) {
                *csv_ << "," << p.branch_mag[static_cast<size_t>(k)]
                      << "," << p.branch_phase_deg[static_cast<size_t>(k)];
            }
            *csv_ << "\n";
        }
        if (retain(r.index)) result.points.push_back(std::move(p));
        if (progress_ && ((r.index + 1) % 100 == 0 || r.index + 1 == nfreq)) {
            progress_(r.index + 1, nfreq, f);
        }
    };

    // ---- solve one block of frequencies ---------------------------------
    // Work is split as groups of `nslots` consecutive frequencies; thread t
    // owns slot t and therefore ALWAYS handles the frequencies congruent to
    // t (mod nslots). That makes the slot assignment deterministic (no
    // dependence on the OpenMP schedule) and race-free, and it is what lets
    // the results be collected back in frequency order afterwards.
    auto solveBlock = [&](int base_index, int items, std::vector<Row>& rows) {
        if (items <= 0) return;
        const int groups = (items + nslots - 1) / nslots;
        long long flops_block = 0;

        if (nslots == 1) {
            ComplexSparseLU& l = lus[0];
            for (int k = 0; k < items; ++k) {
                const size_t idx = static_cast<size_t>(base_index + k);
                const bool good = l.factorize(kPi2 * grid[idx],
                                              g_vals.data(), c_vals.data());
                Row& r = rows[static_cast<size_t>(k)];
                r.index = idx;
                r.valid = good;
                r.re.assign(static_cast<size_t>(nn), 0.0);
                r.im.assign(static_cast<size_t>(nn), 0.0);
                if (!good) continue;
                std::copy(b_re.begin(), b_re.end(), r.re.begin());
                std::copy(b_im.begin(), b_im.end(), r.im.begin());
                l.solve(r.re, r.im);
            }
            flops_block = l.lastFactorFlops() * items;
            if (l.perturbedPivot()) result.stats.pivot_perturbed = true;
        } else {
            std::vector<long long> slot_flops(static_cast<size_t>(nslots), 0);
            std::vector<char> slot_perturbed(static_cast<size_t>(nslots), 0);
            // One thread per slot, each walking the frequencies congruent to
            // its slot: every frequency in [0, items) is covered exactly once,
            // the assignment does not depend on the schedule, and no two
            // threads ever touch the same factor or the same row.
            #pragma omp parallel num_threads(nslots)
            {
                const int t = omp_get_thread_num();
                if (t < nslots) {
                const size_t slot = static_cast<size_t>(t);
                ComplexSparseLU& l = lus[slot];
                for (int g = 0; g < groups; ++g) {
                    const int k = g * nslots + t;
                    if (k >= items) break;
                    const size_t idx = static_cast<size_t>(base_index + k);
                    const bool good = l.factorize(kPi2 * grid[idx],
                                                  g_vals.data(), c_vals.data());
                    slot_flops[slot] += l.lastFactorFlops();
                    if (l.perturbedPivot()) slot_perturbed[slot] = 1;
                    Row& r = rows[static_cast<size_t>(k)];
                    r.index = idx;
                    r.valid = good;
                    r.re.assign(static_cast<size_t>(nn), 0.0);
                    r.im.assign(static_cast<size_t>(nn), 0.0);
                    if (!good) continue;
                    std::copy(b_re.begin(), b_re.end(), r.re.begin());
                    std::copy(b_im.begin(), b_im.end(), r.im.begin());
                    l.solve(r.re, r.im);
                }
                }
            }
            for (int s = 0; s < nslots; ++s) {
                flops_block += slot_flops[static_cast<size_t>(s)];
                if (slot_perturbed[static_cast<size_t>(s)]) {
                    result.stats.pivot_perturbed = true;
                }
            }
        }
        flops_total += flops_block;
        for (const Row& r : rows) {
            if (!r.valid) ++failures_total;
        }
    };

    if (!distribute) {
        // ---------------- local: solve and report block by block ---------
        for (int b_lo = lo_i; b_lo < hi_i; b_lo += static_cast<int>(block)) {
            const int b_hi = std::min(hi_i, b_lo + static_cast<int>(block));
            std::vector<Row> rows(static_cast<size_t>(b_hi - b_lo));
            solveBlock(b_lo, b_hi - b_lo, rows);
            for (const Row& r : rows) report(r);
        }
    } else {
        // ---------------- distributed: exchange in batches ----------------
        // Ranks own contiguous frequency blocks, so a global index maps to a
        // local offset arithmetically (O(1)). Rows travel in chunks of
        // kExchangeChunk_ frequencies: a handful of collectives for the whole
        // sweep instead of one per frequency.
        std::vector<double> local;
        std::vector<double> gathered;
        for (int b_lo = lo_i; b_lo < hi_i; b_lo += static_cast<int>(block)) {
            const int b_hi = std::min(hi_i, b_lo + static_cast<int>(block));
            const int items = b_hi - b_lo;
            std::vector<Row> rows(static_cast<size_t>(items));
            solveBlock(b_lo, items, rows);

            // The block [b_lo, b_hi) is partitioned among the ranks, so every
            // frequency appears exactly once across the gathered chunks. They
            // arrive GROUPED BY RANK (and each rank walks its own chunk
            // offsets), which is not frequency order: collect the whole block
            // and sort it back by global index before reporting.
            std::vector<Row> block_rows;
            block_rows.reserve(static_cast<size_t>(items));
            for (int c_lo = 0; c_lo < items; c_lo += kExchangeChunk_) {
                const int c_hi = std::min(items, c_lo + kExchangeChunk_);
                local.clear();
                local.push_back(static_cast<double>(c_hi - c_lo));
                for (int k = c_lo; k < c_hi; ++k) {
                    const Row& r = rows[static_cast<size_t>(k)];
                    local.push_back(static_cast<double>(r.index));
                    local.push_back(r.valid ? 1.0 : 0.0);
                    local.insert(local.end(), r.re.begin(), r.re.end());
                    local.insert(local.end(), r.im.begin(), r.im.end());
                }
                mpi.gather(local, gathered);
                if (rank != 0) continue;

                size_t off = 0;
                while (off < gathered.size()) {
                    const size_t n = static_cast<size_t>(gathered[off++]);
                    for (size_t k = 0; k < n; ++k) {
                        Row r;
                        r.index = static_cast<size_t>(gathered[off++]);
                        r.valid = gathered[off++] != 0.0;
                        const size_t take = std::min<size_t>(
                            static_cast<size_t>(nn), gathered.size() - off);
                        r.re.assign(gathered.begin() + static_cast<long>(off),
                                    gathered.begin() + static_cast<long>(off + take));
                        off += take;
                        const size_t take2 = std::min<size_t>(
                            static_cast<size_t>(nn), gathered.size() - off);
                        r.im.assign(gathered.begin() + static_cast<long>(off),
                                    gathered.begin() + static_cast<long>(off + take2));
                        off += take2;
                        block_rows.push_back(r);
                    }
                }
            }
            std::sort(block_rows.begin(), block_rows.end(),
                      [](const Row& a, const Row& b) {
                          return a.index < b.index;
                      });
            for (const Row& r : block_rows) {
                report(r);
            }
        }
    }

    // ---- statistics (every rank reports the same numbers) ----------------
    long long total_flops = flops_total;
    int total_failures = failures_total;
    bool perturbed_any = result.stats.pivot_perturbed;
    if (distribute) {
        std::vector<double> stats_pack;
        stats_pack.push_back(static_cast<double>(total_flops));
        stats_pack.push_back(static_cast<double>(total_failures));
        stats_pack.push_back(perturbed_any ? 1.0 : 0.0);
        mpi.allReduceSumInPlace(stats_pack);
        total_flops = static_cast<long long>(stats_pack[0]);
        total_failures = static_cast<int>(stats_pack[1]);
        perturbed_any = stats_pack[2] != 0.0;
    }

    result.stats.frequencies = static_cast<long long>(nfreq);
    result.stats.factorizations = static_cast<int>(nfreq);
    result.stats.factor_flops = total_flops;
    result.stats.failed_frequencies = total_failures;
    result.stats.pivot_perturbed = perturbed_any;
    result.stats.ranks_used = distribute ? world : 1;
    result.stats.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0)
            .count();

    result.success = (total_failures == 0);
    result.message = std::to_string(nfreq) + " frequencies ("
                   + result.variation + "), "
                   + std::to_string(nn) + " unknowns";
    if (total_failures > 0) {
        result.message += ", " + std::to_string(total_failures)
                        + " frequency point(s) could not be factorized";
    }
    if (distribute) {
        result.message += " [frequencies on " + std::to_string(world) + " ranks]";
    }
    if (perturbed_any) {
        result.message += " [pivot perturbation applied]";
    }
    if (result.stats.peak_node >= 0) {
        result.message += ", peak of |" +
                          (result.stats.peak_node < static_cast<int>(csv_labels_.size())
                               ? csv_labels_[static_cast<size_t>(result.stats.peak_node)]
                               : std::string("node?")) +
                          "| = " + std::to_string(result.stats.peak_mag) +
                          " (" + std::to_string(result.stats.peak_db) + " dB) at " +
                          std::to_string(result.stats.peak_frequency) + " Hz";
    }
    return result;
}

} // namespace e_XSpice