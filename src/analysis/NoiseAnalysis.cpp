/**
 * @author Gabriel Henrique Silva
 * @file NoiseAnalysis.cpp
 * @brief Noise analysis: one adjoint solve per frequency gives the transfer of
 *        every noise generator to the output.
 */

#include "NoiseAnalysis.h"

#include "../core/DaeSolution.h"
#include "../core/Device.h"
#include "../core/Node.h"
#include "../parallel/MPIManager.h"
#include "../solver/ComplexSparseLU.h"
#include "SmallSignalLinearization.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace e_XSpice {

namespace {

constexpr double kPi2 = 6.28318530717958647692;   // 2*pi
constexpr double kMinLog  = 1e-38;   // SPICE N_MINLOG: smallest loggable density
constexpr double kMinGain = 1e-20;   // SPICE N_MINGAIN: |gain|^2 floor
constexpr double kIntFlatThresh = 1e-10;   // SPICE N_INTFTHRESH
constexpr double kIntUseLog     = 1e-10;   // SPICE N_INTUSELOG
constexpr double kDegPerRad = 57.29577951308232087680;

/** Case-insensitive comparison (SPICE netlists are not case sensitive). */
bool iequals(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    return true;
}

/**
 * @brief SPICE3's Nintegrate: integral of a·f^k between two sweep points,
 *        fitted in log-log space. Same routine as Xyce's noiseIntegral().
 *
 * The exponent is read off the two densities, so the whole noise spectrum is
 * integrated with a piecewise power law — the same approximation every SPICE
 * uses for "total noise", and the only reason the answer depends slightly on
 * the sweep density.
 */
double noiseIntegral(double dens, double lnDens, double lnLastDens,
                     double delLnFreq, double delFreq, double lnFreq,
                     double lnLastFreq) {
    if (!(delLnFreq > 0.0) || !(dens >= 0.0)) return 0.0;
    const double exponent = (lnDens - lnLastDens) / delLnFreq;
    if (std::fabs(exponent) < kIntFlatThresh) {
        return dens * delFreq;                      // flat density
    }
    const double a = std::exp(lnDens - exponent * lnFreq);
    const double e = exponent + 1.0;
    double integral;
    if (std::fabs(e) < kIntUseLog) {
        integral = a * (lnFreq - lnLastFreq);        // 1/f
    } else {
        integral = a * ((std::exp(e * lnFreq) -
                         std::exp(e * lnLastFreq)) / e);
    }
    if (!std::isfinite(integral)) {
        // A pathologically steep segment would overflow the power-law form;
        // the trapezoid in linear space is the sane fallback and is what a
        // fine sweep would converge to anyway.
        integral = 0.5 * (dens + std::exp(lnLastDens)) * delFreq;
    }
    return integral;
}

} // namespace

// ===============================================================
// Configuration
// ===============================================================

NoiseAnalysis::NoiseAnalysis(Circuit& circuit) : circuit_(circuit) {}
NoiseAnalysis::~NoiseAnalysis() = default;

void NoiseAnalysis::setOutputNode(const std::string& node) {
    output_node_ = node;
    reference_node_ = "0";
}

void NoiseAnalysis::setOutputNode(const std::string& node,
                                  const std::string& ref) {
    output_node_ = node;
    reference_node_ = ref;
}

std::vector<double> NoiseAnalysis::buildGrid(ACVariation variation, int points,
                                             double fstart, double fstop) {
    // Same sweep rule as .AC (and .DC): LIN spreads `points` values between the
    // endpoints, DEC/OCT use `points` per decade/octave and include both ends.
    return ACAnalysis::buildGrid(variation, points, fstart, fstop);
}

// ===============================================================
// Run
// ===============================================================

NoiseResult NoiseAnalysis::run(bool use_mpi) {
    NoiseResult result;
    result.variation = (variation_ == ACVariation::LIN)  ? "LIN"
                     : (variation_ == ACVariation::DEC)  ? "DEC" : "OCT";
    result.output_node = output_node_;
    result.reference_node = reference_node_;
    result.input_source = input_source_;

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
    if (output_node_.empty()) {
        result.message = "no output node given (use V(out) in the .NOISE line)";
        return result;
    }

    MPIManager& mpi = MPIManager::getInstance();
    const int rank = mpi.getRank();
    const int world = mpi.getSize();

    // Frequencies are distributed only when the exchange can pay off.
    const bool distribute = point_distribution_ && use_mpi && world > 1 &&
                            nfreq >= static_cast<size_t>(2) *
                                         static_cast<size_t>(world);

    if (verbose_ && rank == 0) {
        std::cout << "[NoiseAnalysis] " << result.variation << " sweep "
                  << nfreq << " frequencies, " << fstart_ << " Hz .. "
                  << fstop_ << " Hz" << std::endl;
    }

    const auto t0 = std::chrono::high_resolution_clock::now();

    // ---- operating point + linearization (shared with the AC analysis) ----
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
    const int n_full = lin.fullSize();
    const int nn = lin.size();

    if (verbose_ && rank == 0) {
        std::cout << "[NoiseAnalysis] operating point in "
                  << result.stats.newton_iterations << " Newton iterations"
                  << std::endl;
    }

    // ---- resolve the output node(s) ---------------------------------------
    // Returns the DAE row of a node, or -1 for the ground. `found` tells
    // "is ground" (a legal reference) from "does not exist" (an error).
    const auto& node_map = lin.system().nodeMap();
    const std::vector<int>& node_map_for_ctx = node_map;
    auto resolveNode = [&](const std::string& name, bool& found) -> int {
        found = false;
        if (name.empty() || name == "0") return -1;     // ground
        for (const auto& n : circuit_.getNodes()) {
            if (n->isGround()) continue;
            if (!iequals(n->getName(), name)) continue;
            const int id = n->getId();
            if (id < 0 || id >= static_cast<int>(node_map.size())) return -1;
            found = true;
            return node_map[static_cast<size_t>(id)];
        }
        return -1;
    };

    bool out_found = false, ref_found = false;
    const int out_dae = resolveNode(output_node_, out_found);
    const int ref_dae = resolveNode(reference_node_, ref_found);
    if (!out_found) {
        result.success = false;
        result.message = "noise output node '" + output_node_ +
                         "' not found in the circuit";
        if (rank == 0) std::cerr << result.message << std::endl;
        return result;
    }
    if (!ref_found && reference_node_ != "0" && !reference_node_.empty()) {
        result.success = false;
        result.message = "noise reference node '" + reference_node_ +
                         "' not found in the circuit";
        if (rank == 0) std::cerr << result.message << std::endl;
        return result;
    }

    // ---- noise generators at the operating point --------------------------
    // Refresh the device output caches first: in the distributed layout each
    // rank evaluates only its own devices, so a device owned by another rank
    // would report a stale current — and every rank must see the SAME list.
    for (const auto& d : circuit_.getDevices()) {
        d->refreshOutputs(DaeSolution(lin.system().solution(), circuit_,
                                      node_map_for_ctx),
                          0.0);
    }

    std::vector<NoiseSource> raw;
    for (const auto& d : circuit_.getDevices()) {
        raw.clear();
        NoiseContext ctx;
        ctx.temperature = circuit_.temperature();
        ctx.op = nullptr;   // devices report from their own converged caches
        d->loadNoiseSources(ctx, raw);
        for (const auto& s : raw) {
            NoiseSourceInfo info;
            info.device = d->getName();
            info.label = s.label;
            info.pos = (s.pos >= 0) ? lin.acRow(s.pos) : -1;
            info.neg = (s.neg >= 0) ? lin.acRow(s.neg) : -1;
            info.psd_a2 = s.psd_a2;
            info.f_exponent = s.f_exponent;
            // A generator with both terminals grounded has no transfer path.
            if (info.pos < 0 && info.neg < 0) continue;
            if (!(info.psd_a2 > 0.0) && !(info.psd_a2 < 0.0)) continue;
            result.sources.push_back(info);
        }
    }
    const size_t nsrc = result.sources.size();
    result.stats.sources = static_cast<int>(nsrc);

    if (nsrc == 0 && verbose_ && rank == 0) {
        std::cout << "[NoiseAnalysis] no device reports a noise generator"
                  << std::endl;
    }

    // Source list as flat arrays (SoA): the inner loop runs once per frequency
    // over all generators, and flat arrays vectorize where a vector of structs
    // would not.
    std::vector<int> sp_pos(nsrc, -1), sp_neg(nsrc, -1);
    std::vector<double> sp_psd(nsrc, 0.0), sp_fexp(nsrc, 0.0);
    std::vector<unsigned char> sp_white(nsrc, 1);
    for (size_t i = 0; i < nsrc; ++i) {
        sp_pos[i] = result.sources[i].pos;
        sp_neg[i] = result.sources[i].neg;
        sp_psd[i] = result.sources[i].psd_a2;
        sp_fexp[i] = result.sources[i].f_exponent;
        sp_white[i] = (result.sources[i].f_exponent == 0.0) ? 1 : 0;
    }

    // ---- input source (the gain the noise is referred to) -----------------
    // gain = V(out)/V(src), with the input source excited by a UNIT phasor and
    // every other source zeroed. The right-hand side is frequency independent
    // (a source phasor is), so it is assembled once; per frequency this costs
    // ONE forward solve on top of the adjoint solve, and both reuse the SAME
    // factorization — no second factorization, no second matrix build.
    //
    // The MNA matrix is NOT symmetric in general: a source terminal sitting on
    // the ground leaves a +-1 in its branch row with no matching KCL entry
    // once the ground column is eliminated, so the gain cannot be read off the
    // adjoint vector by reciprocity (on a two-resistor divider the reciprocal
    // answer is off by the source's terminal resistance). What the unit
    // excitation DOES buy is independence from the AC magnitude in the
    // netlist: `.NOISE V(out) V1 ...` works even when V1 has no `AC` spec.
    std::vector<double> gain_rhs_re(static_cast<size_t>(nn), 0.0);
    std::vector<double> gain_rhs_im(static_cast<size_t>(nn), 0.0);
    bool have_gain = false;
    {
        std::string want = input_source_;
        Device* src = nullptr;
        if (want.empty()) {
            for (const auto& d : circuit_.getDevices()) {
                if (d->getType() == DeviceType::VOLTAGE_SOURCE) {
                    want = d->getName();
                    src = d.get();
                    break;
                }
            }
        } else {
            for (const auto& d : circuit_.getDevices()) {
                if (iequals(d->getName(), want)) { src = d.get(); break; }
            }
            if (!src) {
                result.success = false;
                result.message = "input source '" + want +
                                 "' not found in the circuit";
                if (rank == 0) std::cerr << result.message << std::endl;
                return result;
            }
        }
        if (src) {
            const DeviceType t = src->getType();
            if (t == DeviceType::VOLTAGE_SOURCE ||
                t == DeviceType::CURRENT_SOURCE) {
                result.input_source = want;
                std::vector<double> full_re(static_cast<size_t>(n_full), 0.0);
                std::vector<double> full_im(static_cast<size_t>(n_full), 0.0);
                // Unit phasor (1∠0) on the selected source, built HERE rather
                // than through loadACExcitation: the noise analysis must not
                // depend on what `AC` magnitude the netlist happens to carry
                // (SPICE/ngspice/Xyce all drive the input with a unit phasor
                // regardless). Sign conventions mirror Device::loadB.
                if (t == DeviceType::VOLTAGE_SOURCE) {
                    const int b = src->getBranchIndex();
                    if (b < 0) {
                        result.success = false;
                        result.message = "input source '" + want +
                                         "' has no branch unknown";
                        if (rank == 0) std::cerr << result.message << std::endl;
                        return result;
                    }
                    full_re[static_cast<size_t>(b)] += 1.0;
                } else {
                    const auto& nodes = src->getNodes();
                    for (size_t k = 0; k < nodes.size() && k < 2; ++k) {
                        const Node* nd = nodes[k];
                        if (!nd) continue;
                        const int id = nd->getId();
                        if (id < 0 ||
                            id >= static_cast<int>(node_map.size())) {
                            continue;
                        }
                        const int row = node_map[static_cast<size_t>(id)];
                        if (row < 0) continue;       // ground
                        // Current source: B[pos] -= I, B[neg] += I
                        full_re[static_cast<size_t>(row)] += (k == 0) ? -1.0
                                                                      : 1.0;
                    }
                }
                gain_rhs_re = lin.permuteRhs(full_re);
                gain_rhs_im = lin.permuteRhs(full_im);
                have_gain = true;
            } else {
                result.success = false;
                result.message = "input source '" + want +
                                 "' is not an independent source";
                if (rank == 0) std::cerr << result.message << std::endl;
                return result;
            }
        } else {
            result.success = false;
            result.message = "no independent source in the circuit: the "
                             "input-referred noise is undefined";
            if (rank == 0) std::cerr << result.message << std::endl;
            return result;
        }
    }

    // ---- adjoint right-hand side (fixed: e_out) ---------------------------
    // Assembled once, in DAE space, then mapped through the row permutation.
    std::vector<double> rhs_full(static_cast<size_t>(n_full), 0.0);
    if (out_dae >= 0) rhs_full[static_cast<size_t>(out_dae)] += 1.0;
    if (ref_dae >= 0) rhs_full[static_cast<size_t>(ref_dae)] -= 1.0;
    const std::vector<double> rhs = lin.permuteRhs(rhs_full);
    std::vector<double> zero_im(static_cast<size_t>(nn), 0.0);
    const int out_ac_row = lin.acRow(out_dae);
    const int ref_ac_row = lin.acRow(ref_dae);

    if (std::getenv("EXSPICE_DEBUG_NOISE")) {
        std::fprintf(stderr, "[NOISE] nn=%d out_row=%d ref_row=%d sources=%zu\n",
                     nn, out_ac_row, ref_ac_row, nsrc);
    }

    // ---- retention --------------------------------------------------------
    size_t stride = 1;
    if (max_retained_ > 0 && nfreq > max_retained_) {
        stride = (nfreq + max_retained_ - 1) / max_retained_;
    }
    const auto retain = [this, stride, nfreq](size_t index) -> bool {
        return max_retained_ == 0 || (index % stride == 0)
               || (index + 1 == nfreq);
    };

    if (csv_) csv_saved_precision_ = csv_->precision(17);

    auto writeHeader = [&] {
        if (!csv_ || csv_header_written_) return;
        *csv_ << "frequency,onoise_psd_v2_per_hz,inoise_psd_v2_per_hz,"
                 "gain_mag,gain_phase_deg";
        if (source_detail_) {
            for (const auto& s : result.sources) {
                *csv_ << "," << s.column() << "_v2_per_hz";
            }
        }
        *csv_ << "\n";
        csv_header_written_ = true;
    };

    // ---- integral accumulators -------------------------------------------
    // SPICE integrates in log space, so the previous point's densities are
    // kept (not the previous point's linear values).
    std::vector<double> source_noise_v2(nsrc, 0.0);
    double out_noise_v2 = 0.0, in_noise_v2 = 0.0;
    std::vector<double> last_ln_src(nsrc, 0.0);
    double last_ln_out = 0.0, last_ln_in = 0.0, last_f = 0.0;
    bool have_previous = false;

    // Per-thread factor slots (see ACAnalysis): frequencies are independent,
    // so each thread owns one factorization workspace and its own buffers.
    int nthreads = 1;
    #ifdef _OPENMP
    nthreads = std::max(1, omp_get_max_threads());
    #endif
    int lo_i = 0, hi_i = static_cast<int>(nfreq);
    if (distribute) mpi.blockRangeForRank(static_cast<int>(nfreq), rank, lo_i, hi_i);
    const int n_local = hi_i - lo_i;
    const bool parallel = (nthreads > 1) && (n_local > 1);
    const int nslots = parallel ? std::min(nthreads, n_local) : 1;

    std::vector<ComplexSparseLU> lus(static_cast<size_t>(nslots));
    // One solution buffer per slot: the adjoint vector is scratch, not a
    // result, so it never leaves the inner loop.
    std::vector<std::vector<double>> scratch_re(static_cast<size_t>(nslots));
    std::vector<std::vector<double>> scratch_im(static_cast<size_t>(nslots));
    std::vector<std::vector<double>> gain_re_(static_cast<size_t>(nslots));
    std::vector<std::vector<double>> gain_im_(static_cast<size_t>(nslots));
    for (int s = 0; s < nslots; ++s) {
        // The symbolic factorization runs once per slot; every frequency then
        // only refills values on this static pattern (the PERMUTED one — it
        // must be the same pattern the value arrays follow).
        lus[static_cast<size_t>(s)].setPattern(lin.rowPointers(),
                                               lin.colIndices());
        scratch_re[static_cast<size_t>(s)].assign(static_cast<size_t>(nn), 0.0);
        scratch_im[static_cast<size_t>(s)].assign(static_cast<size_t>(nn), 0.0);
        gain_re_[static_cast<size_t>(s)].assign(static_cast<size_t>(nn), 0.0);
        gain_im_[static_cast<size_t>(s)].assign(static_cast<size_t>(nn), 0.0);
    }

    // Frequencies are processed in BLOCKS so RAM stays bounded.
    const size_t payload_cap = static_cast<size_t>(1) << 20;
    // The per-generator PSDs are ALWAYS computed, not only when the CSV asks
    // for the columns: the per-generator ranking is the practical output of a
    // noise analysis ("what do I fix first") and it must not depend on a
    // printing option. `pts_per_summary` controls the CSV columns only.
    const size_t per_freq = 5 + nsrc;
    const size_t block = std::max<size_t>(
        1, std::min<size_t>(4096,
                             payload_cap / std::max<size_t>(1, per_freq)));

    // Per-frequency payload: [valid, onoise, inoise, gain_re, gain_im]
    constexpr size_t kScalars = 5;
    std::vector<double> scalar_data;
    std::vector<double> psd_data;
    std::vector<unsigned char> row_valid;
    long long flops_total = 0;
    int failures_total = 0;

    // ---- report (always in frequency order: the integral depends on it) ---
    auto report = [&](size_t index, bool valid, const double* scalars,
                      const double* psd) {
        if (!valid) {
            ++failures_total;
            return;
        }
        const double f = grid[index];
        const double onoise = scalars[1];
        const double inoise = scalars[2];
        const double gre = scalars[3], gim = scalars[4];

        writeHeader();
        if (csv_) {
            *csv_ << f << "," << onoise << "," << inoise << ","
                  << std::hypot(gre, gim) << "," << std::atan2(gim, gre) * kDegPerRad;
            if (source_detail_) {
                for (size_t i = 0; i < nsrc; ++i) *csv_ << "," << psd[i];
            }
            *csv_ << "\n";
        }

        if (retain(index)) {
            NoiseDataPoint p;
            p.frequency = f;
            p.onoise_psd = onoise;
            p.inoise_psd = inoise;
            p.gain_mag = std::hypot(gre, gim);
            p.gain_phase_deg = std::atan2(gim, gre) * kDegPerRad;
            p.onoise_amp = std::sqrt(std::max(onoise, 0.0));
            p.inoise_amp = std::sqrt(std::max(inoise, 0.0));
            p.valid = true;
            if (psd) {
                // Kept only when the CSV asked for the columns: the retained
                // sample is a decimated log aid, and n_sources x points would
                // defeat the RAM bound the block loop exists for. The noise
                // BUDGET (per-generator integrals) is accumulated below for
                // every generator regardless.
                p.source_psd.assign(psd, psd + nsrc);
            }
            result.points.push_back(std::move(p));
        }

        if (onoise > result.stats.peak_output_psd) {
            result.stats.peak_output_psd = onoise;
            result.stats.peak_output_freq = f;
        }

        // ---- integration (SPICE Nintegrate, log-log power law) ----------
        const double lnF = std::log(f);
        if (have_previous) {
            const double lnLastF = std::log(last_f);
            const double delLnF = lnF - lnLastF;
            const double delF = f - last_f;
            const double lnOut = std::log(std::max(onoise, kMinLog));
            const double lnIn = std::log(std::max(inoise, kMinLog));
            out_noise_v2 += noiseIntegral(onoise, lnOut, last_ln_out, delLnF,
                                          delF, lnF, lnLastF);
            in_noise_v2 += noiseIntegral(inoise, lnIn, last_ln_in, delLnF,
                                         delF, lnF, lnLastF);
            for (size_t i = 0; i < nsrc; ++i) {
                if (!psd) break;
                const double d = psd[i];
                source_noise_v2[i] += noiseIntegral(
                    d, std::log(std::max(d, kMinLog)), last_ln_src[i],
                    delLnF, delF, lnF, lnLastF);
            }
        }
        last_ln_out = std::log(std::max(onoise, kMinLog));
        last_ln_in = std::log(std::max(inoise, kMinLog));
        for (size_t i = 0; i < nsrc; ++i) {
            last_ln_src[i] = std::log(std::max(psd ? psd[i] : 0.0, kMinLog));
        }
        last_f = f;
        have_previous = true;

        if (progress_ && ((index + 1) % 100 == 0 || index + 1 == nfreq)) {
            progress_(index + 1, nfreq, f);
        }
    };

    // ---- solve one block of frequencies ------------------------------------
    auto solveBlock = [&](int base_index, int items, double* scalars,
                          double* psd) {
        if (items <= 0) return;
        const int groups = (items + nslots - 1) / nslots;
        long long flops_block = 0;
        const double* g_ptr = lin.gValues().data();
        const double* c_ptr = lin.cValues().data();

        auto solveOne = [&](ComplexSparseLU& lu, int slot, int k) {
            const size_t idx = static_cast<size_t>(base_index + k);
            const double f = grid[idx];
            double* sc = scalars + static_cast<size_t>(k) * kScalars;
            sc[0] = 0.0;
            if (!lu.factorize(kPi2 * f, g_ptr, c_ptr)) return;
            std::vector<double>& xr = scratch_re[static_cast<size_t>(slot)];
            std::vector<double>& xi = scratch_im[static_cast<size_t>(slot)];
            std::copy(rhs.begin(), rhs.end(), xr.begin());
            std::copy(zero_im.begin(), zero_im.end(), xi.begin());
            lu.solveTranspose(xr, xi);

            double total = 0.0;
            double* p = psd ? psd + static_cast<size_t>(k) * nsrc : nullptr;
            for (size_t i = 0; i < nsrc; ++i) {
                const int a = sp_pos[i], b = sp_neg[i];
                const double hr = (a >= 0 ? xr[static_cast<size_t>(a)] : 0.0)
                                - (b >= 0 ? xr[static_cast<size_t>(b)] : 0.0);
                const double hi = (a >= 0 ? xi[static_cast<size_t>(a)] : 0.0)
                                - (b >= 0 ? xi[static_cast<size_t>(b)] : 0.0);
                const double s = sp_white[i]
                    ? sp_psd[i]
                    : sp_psd[i] * std::pow(f, -sp_fexp[i]);
                const double contrib = (hr * hr + hi * hi) * s;
                if (p) p[i] = contrib;
                total += contrib;
            }
            double gre = 0.0, gim = 0.0;
            if (have_gain) {
                // Gain from the SAME factorization: one forward solve with the
                // unit phasor of the input source.
                std::vector<double>& ar = gain_re_[static_cast<size_t>(slot)];
                std::vector<double>& ai = gain_im_[static_cast<size_t>(slot)];
                std::copy(gain_rhs_re.begin(), gain_rhs_re.end(), ar.begin());
                std::copy(gain_rhs_im.begin(), gain_rhs_im.end(), ai.begin());
                lu.solve(ar, ai);
                const int ro = out_ac_row, rr = ref_ac_row;
                gre = (ro >= 0 ? ar[static_cast<size_t>(ro)] : 0.0)
                    - (rr >= 0 ? ar[static_cast<size_t>(rr)] : 0.0);
                gim = (ro >= 0 ? ai[static_cast<size_t>(ro)] : 0.0)
                    - (rr >= 0 ? ai[static_cast<size_t>(rr)] : 0.0);
            }
            const double gain2 = gre * gre + gim * gim;
            sc[1] = total;
            sc[2] = have_gain ? total / std::max(gain2, kMinGain) : 0.0;
            sc[3] = gre;
            sc[4] = gim;
            sc[0] = 1.0;
        };

        if (nslots == 1) {
            for (int k = 0; k < items; ++k) solveOne(lus[0], 0, k);
            flops_block = lus[0].lastFactorFlops() * items;
            if (lus[0].perturbedPivot()) result.stats.pivot_perturbed = true;
        } else {
            // Each slot accumulates the flops of ITS OWN factorizations: the
            // LU only remembers the last one, so reading it after each solve
            // is what makes the total match the serial path.
            std::vector<long long> slot_flops(static_cast<size_t>(nslots), 0);
            #pragma omp parallel num_threads(nslots)
            {
                const int t = omp_get_thread_num();
                if (t < nslots) {
                    const size_t slot = static_cast<size_t>(t);
                    for (int g = 0; g < groups; ++g) {
                        const int k = g * nslots + t;
                        if (k >= items) break;
                        solveOne(lus[slot], static_cast<int>(slot), k);
                        slot_flops[slot] += lus[slot].lastFactorFlops();
                    }
                }
            }
            for (int s = 0; s < nslots; ++s) {
                flops_block += slot_flops[static_cast<size_t>(s)];
                if (lus[static_cast<size_t>(s)].perturbedPivot()) {
                    result.stats.pivot_perturbed = true;
                }
            }
        }
        flops_total += flops_block;
    };

    if (!distribute) {
        for (int b_lo = lo_i; b_lo < hi_i; b_lo += static_cast<int>(block)) {
            const int b_hi = std::min(hi_i, b_lo + static_cast<int>(block));
            const int items = b_hi - b_lo;
            scalar_data.assign(static_cast<size_t>(items) * kScalars, 0.0);
            psd_data.assign(static_cast<size_t>(items) * nsrc, 0.0);
            solveBlock(b_lo, items, scalar_data.data(), psd_data.data());
            for (int k = 0; k < items; ++k) {
                const double* sc = scalar_data.data() +
                                   static_cast<size_t>(k) * kScalars;
                const double* p = psd_data.data() +
                                  static_cast<size_t>(k) * nsrc;
                report(static_cast<size_t>(b_lo + k), sc[0] != 0.0, sc, p);
            }
        }
    } else {
        std::vector<double> local;
        std::vector<double> gathered;
        size_t emitted = 0;
        for (int b_lo = lo_i; b_lo < hi_i; b_lo += static_cast<int>(block)) {
            const int b_hi = std::min(hi_i, b_lo + static_cast<int>(block));
            const int items = b_hi - b_lo;
            scalar_data.assign(static_cast<size_t>(items) * kScalars, 0.0);
            psd_data.assign(static_cast<size_t>(items) * nsrc, 0.0);
            solveBlock(b_lo, items, scalar_data.data(), psd_data.data());

            for (int c_lo = 0; c_lo < items; c_lo += kExchangeChunk_) {
                const int c_hi = std::min(items, c_lo + kExchangeChunk_);
                local.clear();
                local.push_back(static_cast<double>(c_hi - c_lo));
                for (int k = c_lo; k < c_hi; ++k) {
                    const double* sc = scalar_data.data() +
                                       static_cast<size_t>(k) * kScalars;
                    local.push_back(static_cast<double>(b_lo + k));
                    local.push_back(sc[0]);
                    for (size_t s = 1; s < kScalars; ++s) local.push_back(sc[s]);
                    // The per-generator PSDs travel with the row even when the
                    // CSV does not print them: the noise budget (per-generator
                    // integrals) is computed on the master, which needs them.
                    const double* p = psd_data.data() +
                                      static_cast<size_t>(k) * nsrc;
                    local.insert(local.end(), p, p + nsrc);
                }
                mpi.gather(local, gathered);
                if (rank != 0) continue;

                size_t off = 0;
                while (off < gathered.size()) {
                    const size_t n = static_cast<size_t>(gathered[off++]);
                    for (size_t k = 0; k < n; ++k) {
                        const size_t index =
                            static_cast<size_t>(gathered[off++]);
                        const double valid = gathered[off++];
                        double sc[kScalars] = {0.0, 0.0, 0.0, 0.0, 0.0};
                        for (size_t s = 1; s < kScalars; ++s) sc[s] = gathered[off++];
                        const double* p = &gathered[off];
                        off += nsrc;
                        if (index != emitted) continue;   // defensive
                        report(index, valid != 0.0, sc, p);
                        ++emitted;
                    }
                }
            }
        }
    }

    // ---- reduce the per-rank partial sums ----------------------------------
    long long total_flops = flops_total;
    int total_failures = failures_total;
    bool perturbed_any = result.stats.pivot_perturbed;
    if (distribute) {
        std::vector<double> pack;
        pack.push_back(static_cast<double>(total_flops));
        pack.push_back(static_cast<double>(total_failures));
        pack.push_back(perturbed_any ? 1.0 : 0.0);
        // The integrals are a sum over the whole sweep, so they reduce like any
        // other accumulator; the ORDER of the sweep does not change their sum
        // (same terms, same order within a rank's contiguous block).
        pack.push_back(out_noise_v2);
        pack.push_back(in_noise_v2);
        pack.push_back(result.stats.peak_output_psd);
        for (double v : source_noise_v2) pack.push_back(v);
        mpi.allReduceSumInPlace(pack);
        total_flops = static_cast<long long>(pack[0]);
        total_failures = static_cast<int>(pack[1]);
        perturbed_any = pack[2] != 0.0;
        out_noise_v2 = pack[3];
        in_noise_v2 = pack[4];
        result.stats.peak_output_psd = pack[5];
        for (size_t i = 0; i < nsrc; ++i) source_noise_v2[i] = pack[6 + i];
        // The peak frequency travels with the peak density: re-scan the rows
        // already emitted (the master has them all) instead of shipping it.
        for (const auto& p : result.points) {
            if (p.onoise_psd == result.stats.peak_output_psd) {
                result.stats.peak_output_freq = p.frequency;
                break;
            }
        }
    }

    result.stats.frequencies = static_cast<long long>(nfreq);
    result.stats.factorizations = static_cast<int>(nfreq);
    result.stats.factor_flops = total_flops;
    result.stats.failed_frequencies = total_failures;
    result.stats.pivot_perturbed = perturbed_any;
    result.stats.ranks_used = distribute ? world : 1;
    result.stats.output_noise_v2 = out_noise_v2;
    result.stats.input_noise_v2 = in_noise_v2;
    result.stats.output_noise_rms = std::sqrt(std::max(out_noise_v2, 0.0));
    result.stats.input_noise_rms = std::sqrt(std::max(in_noise_v2, 0.0));
    result.stats.source_noise_v2 = source_noise_v2;
    result.stats.elapsed_seconds =
        std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();

    result.success = (total_failures == 0);
    result.message = std::to_string(nfreq) + " frequencies ("
                   + result.variation + "), " + std::to_string(nn)
                   + " unknowns, " + std::to_string(nsrc) + " noise source(s)";
    if (total_failures > 0) {
        result.message += ", " + std::to_string(total_failures)
                        + " frequency point(s) could not be factorized";
    }
    if (distribute) {
        result.message += " [frequencies on " + std::to_string(world)
                        + " ranks]";
    }
    if (perturbed_any) {
        result.message += " [pivot perturbation applied]";
    }
    result.message += ", output noise " + std::to_string(result.stats.output_noise_rms)
                    + " V rms, input noise "
                    + std::to_string(result.stats.input_noise_rms) + " V rms";
    return result;
}

} // namespace e_XSpice