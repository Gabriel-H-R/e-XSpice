/**
 * @author Gabriel Henrique Silva
 * @file main.cpp
 * @brief Main entry point of the circuit simulator (file mode only)
 *
 * Usage:
 *   e-XSpice file <netlist.cir> [output.csv]
 *
 * The netlist selects the analysis (.OP, .TRAN, .DC or .AC). Node voltages are
 * reported with their netlist names everywhere (log and CSV): a node
 * called `out` in the netlist appears as `Vout`.
 */

#include "core/Circuit.h"
#include "core/Node.h"
#include "core/NetlistParser.h"
#include "analysis/ACAnalysis.h"
#include "analysis/DCAnalysis.h"
#include "analysis/DCSweep.h"
#include "analysis/NoiseAnalysis.h"
#include "analysis/TFAnalysis.h"
#include "analysis/TransientAnalysis.h"
#include "parallel/MPIManager.h"

#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace e_XSpice;

/**
 * @brief Cap on DC sweep points kept in RAM for the log sample. The CSV
 *        stream keeps every point, so this bounds memory on long sweeps
 *        without losing data on disk.
 */
static constexpr size_t kDCMaxRetainedPoints = 20000;

/**
 * @brief Same idea for the AC sweep: the CSV keeps every frequency, RAM keeps
 *        a decimated sample.
 */
static constexpr size_t kACMaxRetainedPoints = 20000;

/**
 * @brief Builds the `V<name>` labels for every non-ground node, aligned
 *        by position with Circuit::getNodes() — and therefore with
 *        DCResult::node_voltages and TransientDataPoint::node_voltages,
 *        which are collected in that same order.
 */
static std::vector<std::string> nodeLabels(const Circuit& circuit) {
    std::vector<std::string> labels;
    labels.reserve(circuit.getNodes().size());
    for (const auto& n : circuit.getNodes()) {
        const std::string name = n->getName();
        labels.push_back("V" + (name.empty()
                                ? "N" + std::to_string(n->getId())
                                : name));
    }
    return labels;
}

/**
 * @brief Prints one AC log line, magnitude AND decibels.
 *
 * Both are printed on purpose: a frequency response is read in dB (0 dB is
 * the input level, so a roll-off shows up as negative dB and a resonance as a
 * positive bump), while the linear magnitude alone decays "towards zero" and
 * leaves it unclear whether the response ever crossed unity. Phases are in
 * degrees, matching the CSV.
 */
static void printACPoint(const ACDataPoint& p,
                         const std::vector<std::string>& labels,
                         const std::string& prefix = "  ") {
    std::cout << prefix << "f=" << p.frequency << " Hz: ";
    for (size_t j = 0; j < p.mag.size() && j < labels.size(); ++j) {
        std::cout << "|" << labels[j] << "|=" << p.mag[j]
                  << " (" << p.db[j] << " dB) @" << p.phase_deg[j] << "deg ";
    }
    for (size_t k = 0; k < p.branch_mag.size(); ++k) {
        const double db = 20.0 * std::log10(std::max(p.branch_mag[k], 1e-300));
        std::cout << "|I|=" << p.branch_mag[k] << " (" << db << " dB) @"
                  << p.branch_phase_deg[k] << "deg ";
    }
    if (!p.valid) std::cout << "(no solution) ";
    std::cout << std::endl;
}

/**
 * @brief Prints one log line with every node voltage:
 *        `OP: Vout=5V V1=...V ...` (DC) or `t=...s: Vout=...V ...` (TRAN).
 */
static void printVoltageLine(const std::string& prefix,
                             const std::vector<std::string>& labels,
                             const std::vector<double>& voltages) {
    std::cout << "  " << prefix << ": ";
    for (size_t j = 0; j < labels.size() && j < voltages.size(); ++j) {
        if (j > 0) std::cout << " ";
        std::cout << labels[j] << "=" << voltages[j] << "V";
    }
    std::cout << std::endl;
}

/**
 * @brief Simulates circuit from netlist file
 */
bool simulateNetlist(const std::string& filename, bool use_mpi,
                     const std::string& csv_path) {
    MPIManager& mpi = MPIManager::getInstance();

    Circuit circuit("Netlist Circuit");
    NetlistParser parser;
    parser.setVerbose(mpi.isMaster());

    if (mpi.isMaster()) std::cout << "\nParsing netlist: " << filename << std::endl;

    if (!parser.parse(filename, circuit)) {
        if (mpi.isMaster()) std::cerr << "Parse error: " << parser.getErrorMessage() << std::endl;
        return false;
    }

    if (mpi.isMaster()) circuit.printTopology();

    const auto& params = parser.getAnalysisParameters();
    const std::vector<std::string> labels = nodeLabels(circuit);

    if (mpi.isMaster() && !params.option_warnings.empty()) {
        std::cout << "\nWarning: unsupported .OPTIONS entries ignored:";
        for (const auto& w : params.option_warnings) std::cout << " " << w;
        std::cout << std::endl;
    }

    // `.TEMP <celsius>` / `.OPTIONS TEMP=`: re-scale every model card that
    // did not pin its own TEMP= and refresh the devices that precompute
    // temperature scalars. Applied BEFORE any analysis so the first operating
    // point is already at temperature (a .DC TEMP sweep moves it later).
    if (params.temp_given) {
        circuit.setTemperature(params.temp_value + 273.15);
        if (mpi.isMaster()) {
            std::cout << "\nCircuit temperature: " << params.temp_value
                      << " C (" << circuit.temperature() << " K)" << std::endl;
        }
    }

    // A netlist may request both analyses (.OP then .TRAN, in any order):
    // report the operating point first, then run the transient (which
    // re-solves its own internal OP as the initial state).
    const bool do_op = params.want_op || params.type == AnalysisType::DC_OP;
    const bool do_tran = params.want_tran || params.type == AnalysisType::TRANSIENT;
    const bool do_dc = params.want_dc_sweep || params.type == AnalysisType::DC_SWEEP;
    const bool do_ac = params.want_ac || params.type == AnalysisType::AC;
    const bool do_noise = params.want_noise || params.type == AnalysisType::NOISE;
    const bool do_tf = params.want_tf || params.type == AnalysisType::TF;

    if (do_op) {
        DCAnalysis dc(circuit);
        dc.setVerbose(mpi.isMaster());
        dc.setTolerance(1e-10);

        DCResult result = dc.run(use_mpi, false);

        if (!result.success) {
            if (mpi.isMaster()) std::cerr << "DC Analysis failed: " << result.message << std::endl;
            mpi.barrier();
            return false;
        }

        if (mpi.isMaster()) {
            printVoltageLine("OP", labels, result.node_voltages);
            circuit.printSolution();
        }
    }

    if (do_tran) {
        TransientAnalysis tran(circuit);
        tran.setVerbose(mpi.isMaster());
        tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
        tran.setTolerance(1e-9);
        tran.setSaveInterval(10);

        if (mpi.isMaster()) {
            tran.setProgressCallback([](double t, int step, int total) {
                if (step % 1000 == 0) {
                    std::cout << "  Progress: " << (100.0 * step / total) << "% (t=" << t << "s)" << std::endl;
                }
            });
        }

        TransientResult result = tran.run(params.tran_tstart, params.tran_tstop, params.tran_tstep, use_mpi);

        if (!result.success) {
            if (mpi.isMaster()) std::cerr << "Transient Analysis failed: " << result.message << std::endl;
            mpi.barrier();
            return false;
        }

        if (mpi.isMaster()) {
            std::cout << "\nTransient analysis completed!" << std::endl;
            std::cout << "Total simulation time: " << result.total_time << " s" << std::endl;

            int sample_interval = std::max(1, (int)result.data.size() / 10);
            std::cout << "\nSample data points:" << std::endl;
            for (size_t i = 0; i < result.data.size(); i += sample_interval) {
                const auto& point = result.data[i];
                std::cout << "  t=" << point.time << "s: ";
                for (size_t j = 0; j < std::min(size_t(3), point.node_voltages.size()); ++j) {
                    std::cout << labels[j] << "=" << point.node_voltages[j] << "V ";
                }
                std::cout << std::endl;
            }

            // --- Data Export to CSV (netlist node names as headers) ---
            if (!result.data.empty()) {
                std::ofstream csv(csv_path);
                if (csv.is_open()) {
                    csv.precision(17);
                    csv << "time";
                    for (const auto& label : labels) {
                        csv << "," << label;
                    }
                    csv << "\n";

                    for (const auto& point : result.data) {
                        csv << point.time;
                        for (double v : point.node_voltages) {
                            csv << "," << v;
                        }
                        csv << "\n";
                    }
                    csv.close();
                    std::cout << "\n[e-XSpice] Dados salvos em: " << csv_path << "\n";
                } else {
                    std::cerr << "Erro ao criar o arquivo " << csv_path << "\n";
                }
            }
        }
    } else if (!do_op && !do_dc && !do_ac && !do_noise && !do_tf) {
        if (mpi.isMaster()) std::cerr << "No analysis specified in netlist (use .OP, .TRAN, .DC, .AC or .NOISE)" << std::endl;
        return false;
    }

    if (do_dc) {
        auto variationOf = [](const std::string& v) {
            return (v == "DEC") ? DCSweepVariation::DEC
                 : (v == "OCT") ? DCSweepVariation::OCT
                                : DCSweepVariation::LIN;
        };

        DCSweepSpec spec;
        spec.source    = params.dc_source;
        spec.variation = variationOf(params.dc_variation);
        spec.start     = params.dc_start;
        spec.stop      = params.dc_stop;
        spec.step      = params.dc_step;

        DCSweep sweep(circuit);
        sweep.setSweep(spec);
        sweep.setVerbose(mpi.isMaster());
        sweep.setTolerance(1e-10);
        sweep.setContinueOnFail(params.opt_nofail);
        sweep.setMaxRetries(params.opt_max_retries);
        // `.OPTIONS SWEEPPAR`: split the sweep points across ranks instead of
        // running one dependent chain everywhere.
        sweep.setPointDistribution(params.opt_sweep_parallel);
        if (params.opt_max_iterations > 0) {
            sweep.setMaxIterations(params.opt_max_iterations);
        }
        if (params.dc_nested) {
            DCSweepSpec inner;
            inner.source    = params.dc_source2;
            inner.variation = variationOf(params.dc_variation2);
            inner.start     = params.dc_start2;
            inner.stop      = params.dc_stop2;
            inner.step      = params.dc_step2;
            sweep.setInnerSweep(inner);
        }

        // Data Export to CSV, streamed row by row during the sweep: a long
        // sweep no longer holds points x nodes in RAM, and the decimated
        // in-memory sample is only used for the log lines below.
        std::ofstream csv;
        if (mpi.isMaster() && !csv_path.empty()) {
            csv.open(csv_path);
            if (csv.is_open()) {
                csv.precision(17);
                sweep.setCsvStream(&csv, labels);
            } else {
                std::cerr << "Erro ao criar o arquivo " << csv_path << "\n";
            }
        }
        sweep.setMaxRetainedPoints(kDCMaxRetainedPoints);

        if (mpi.isMaster()) {
            sweep.setProgressCallback([](size_t done, size_t total, double v) {
                if (done == total || done % 1000 == 0) {
                    std::cout << "  DC sweep: " << done << "/" << total
                              << " points (last value " << v << ")"
                              << std::endl;
                }
            });
        }

        DCSweepResult result = sweep.run(use_mpi);

        if (mpi.isMaster() && csv.is_open()) {
            csv.close();
            std::cout << "\n[e-XSpice] Dados salvos em: " << csv_path << "\n";
        }

        if (!result.success) {
            if (mpi.isMaster()) std::cerr << "DC Sweep failed: " << result.message << std::endl;
            mpi.barrier();
            return false;
        }

        if (mpi.isMaster()) {
            std::cout << "\nDC sweep completed: " << result.message << std::endl;
            std::cout << "  Newton iterations: " << result.stats.newton_iterations
                      << " | linear solves: " << result.stats.linear_solves
                      << " | predictor points: " << result.stats.predictor_hits
                      << " | source-step retries: " << result.stats.retry_solves
                      << " | ranks: " << result.stats.ranks_used
                      << " | time: " << result.stats.elapsed_seconds << " s"
                      << std::endl;

            int sample_interval = std::max(1, (int)result.points.size() / 10);
            std::cout << "\nSample data points:" << std::endl;
            for (size_t i = 0; i < result.points.size(); i += sample_interval) {
                const auto& point = result.points[i];
                std::cout << "  ";
                if (result.nested) std::cout << result.sweep_source << "=" << point.outer_value << " ";
                // Parenthesized: `<<` binds tighter than `?:`, so without them this would
                // print a boolean instead of the swept variable name.
                std::cout << (result.sweep_source2.empty() ? result.sweep_source
                                                          : result.sweep_source2);
                std::cout << "=" << point.sweep_value << ": ";
                for (size_t j = 0; j < std::min(size_t(3), point.node_voltages.size()); ++j) {
                    std::cout << labels[j] << "=" << point.node_voltages[j] << "V ";
                }
                if (!point.converged) std::cout << "(no solution) ";
                std::cout << std::endl;
            }
        }
    }

    if (do_ac) {
        // Small-signal sweep: ONE operating point + one linearization, then an
        // independent complex solve per frequency (frequencies can be solved
        // in parallel because they do not depend on each other).
        ACAnalysis ac(circuit);
        ac.setVariation((params.ac_variation == "LIN")  ? ACVariation::LIN
                      : (params.ac_variation == "OCT") ? ACVariation::OCT
                                                       : ACVariation::DEC);
        ac.setPoints(params.ac_points);
        ac.setFrequencyRange(params.ac_fstart, params.ac_fstop);
        ac.setTolerance(1e-10);
        ac.setVerbose(mpi.isMaster());
        if (params.opt_max_iterations > 0) {
            ac.setMaxIterations(params.opt_max_iterations);
        }
        // `.OPTIONS ACFREQPAR`: split the FREQUENCIES across ranks. Results
        // stay bit-identical (each rank factorizes the same replicated G, C).
        ac.setPointDistribution(params.opt_ac_freq_parallel);

        std::ofstream csv;
        if (mpi.isMaster() && !csv_path.empty()) {
            csv.open(csv_path);
            if (csv.is_open()) {
                csv.precision(17);
                ac.setCsvStream(&csv, labels);
            } else {
                std::cerr << "Erro ao criar o arquivo " << csv_path << "\n";
            }
        }
        ac.setMaxRetainedPoints(kACMaxRetainedPoints);
        if (mpi.isMaster()) {
            ac.setProgressCallback([](size_t done, size_t total, double f) {
                std::cout << "  AC sweep: " << done << "/" << total
                          << " frequencies (last " << f << " Hz)" << std::endl;
            });
        }

        ACResult result = ac.run(use_mpi);

        if (mpi.isMaster() && csv.is_open()) {
            csv.close();
            std::cout << "\n[e-XSpice] Dados salvos em: " << csv_path << "\n";
        }

        if (!result.success) {
            if (mpi.isMaster()) std::cerr << "AC analysis failed: " << result.message << std::endl;
            mpi.barrier();
            return false;
        }

        if (mpi.isMaster()) {
            std::cout << "\nAC analysis completed: " << result.message << std::endl;
            std::cout << "  OP iterations: " << result.stats.newton_iterations
                      << " | factorizations: " << result.stats.factorizations
                      << " | LU flops: " << result.stats.factor_flops
                      << " | ranks: " << result.stats.ranks_used
                      << " | time: " << result.stats.elapsed_seconds << " s"
                      << std::endl;

            // A frequency response is read in dB (0 dB = the input level, and
            // a roll-off shows as negative dB); the linear magnitude alone
            // decays "towards zero" and hides whether the circuit ever
            // crossed unity. Both are printed, plus the peak of every node so
            // a narrow resonance cannot be missed by the sampling below.
            if (!result.stats.node_peak_mag.empty()) {
                std::cout << "  Peaks (of the sampled grid — a very high-Q "
                             "resonance can be narrower than the frequency step):"
                          << std::endl;
                const size_t shown = std::min<size_t>(4, result.stats.node_peak_mag.size());
                for (size_t i = 0; i < shown; ++i) {
                    const std::string& nm = i < labels.size() ? labels[i] : std::string("?");
                    std::cout << "    |" << nm << "|_max = "
                              << result.stats.node_peak_mag[i] << " ("
                              << result.stats.node_peak_db[i] << " dB) at "
                              << result.stats.node_peak_freq[i] << " Hz" << std::endl;
                }
                if (result.stats.node_peak_mag.size() > shown) {
                    std::cout << "    (" << result.stats.node_peak_mag.size() - shown
                              << " more node(s) in the CSV)" << std::endl;
                }
            }

            // Sample: the retained points, uniformly, PLUS the peak row so a
            // sharp resonance is always visible in the log.
            const size_t stride = std::max<size_t>(
                1, (result.points.size() + 9) / 10);
            std::cout << "\nMagnitude / phase (linear and dB):" << std::endl;
            for (size_t i = 0; i < result.points.size(); i += stride) {
                printACPoint(result.points[i], labels);
            }
            if (result.stats.peak_node >= 0 && !result.points.empty()) {
                // Print the row at the peak frequency reported above. The
                // search is by FREQUENCY, not by magnitude: stats already holds
                // the exact peak over the whole sweep, and the retained sample
                // may be decimated (then this is the nearest row to it).
                const size_t pn = static_cast<size_t>(result.stats.peak_node);
                const double f_peak =
                    pn < result.stats.node_peak_freq.size()
                        ? result.stats.node_peak_freq[pn]
                        : 0.0;
                size_t best = 0;
                for (size_t k = 1; k < result.points.size(); ++k) {
                    if (std::fabs(result.points[k].frequency - f_peak) <
                        std::fabs(result.points[best].frequency - f_peak)) {
                        best = k;
                    }
                }
                if (best % stride != 0) {
                    printACPoint(result.points[best], labels, "  (peak) ");
                }
            }
        }
    }

    if (do_noise) {
        auto variationOf = [](const std::string& v) {
            return (v == "DEC") ? ACVariation::DEC
                 : (v == "OCT") ? ACVariation::OCT
                                : ACVariation::LIN;
        };

        NoiseAnalysis noise(circuit);
        noise.setVariation(variationOf(params.ac_variation));
        noise.setPoints(params.ac_points);
        noise.setFrequencyRange(params.ac_fstart, params.ac_fstop);
        if (params.noise_ref_node == "0" || params.noise_ref_node.empty()) {
            noise.setOutputNode(params.noise_out_node);
        } else {
            noise.setOutputNode(params.noise_out_node, params.noise_ref_node);
        }
        if (params.noise_source_given) noise.setInputSource(params.noise_source);
        // SPICE's pts_per_summary: per-generator contribution breakdown.
        noise.setSourceDetail(params.noise_detail);
        noise.setTolerance(1e-10);
        noise.setVerbose(mpi.isMaster());
        if (params.opt_max_iterations > 0) {
            noise.setMaxIterations(params.opt_max_iterations);
        }
        // `.OPTIONS ACFREQPAR` splits the frequencies across ranks: each rank
        // runs the SAME adjoint solve on its own frequencies, so the answer is
        // bit-identical to a serial run.
        noise.setPointDistribution(params.opt_ac_freq_parallel);

        std::ofstream csv;
        if (mpi.isMaster() && !csv_path.empty()) {
            csv.open(csv_path);
            if (csv.is_open()) {
                csv.precision(17);
                noise.setCsvStream(&csv);
            } else {
                std::cerr << "Erro ao criar o arquivo " << csv_path << "\n";
            }
        }
        noise.setMaxRetainedPoints(kACMaxRetainedPoints);
        if (mpi.isMaster()) {
            noise.setProgressCallback([](size_t done, size_t total, double f) {
                std::cout << "  Noise sweep: " << done << "/" << total
                          << " frequencies (last " << f << " Hz)" << std::endl;
            });
        }

        NoiseResult result = noise.run(use_mpi);

        if (mpi.isMaster() && csv.is_open()) {
            csv.close();
            std::cout << "\n[e-XSpice] Dados salvos em: " << csv_path << "\n";
        }

        if (!result.success) {
            if (mpi.isMaster()) {
                std::cerr << "Noise analysis failed: " << result.message
                          << std::endl;
            }
            mpi.barrier();
            return false;
        }

        if (mpi.isMaster()) {
            std::cout << "\nNoise analysis completed: " << result.message
                      << std::endl;
            std::cout << "  OP iterations: " << result.stats.newton_iterations
                      << " | factorizations: " << result.stats.factorizations
                      << " | LU flops: " << result.stats.factor_flops
                      << " | ranks: " << result.stats.ranks_used
                      << " | time: " << result.stats.elapsed_seconds << " s"
                      << std::endl;

            // The integrated numbers are what a noise budget is quoted in
            // (rms over the analysed band), so they lead.
            std::cout << "\nIntegrated noise over the swept band:" << std::endl;
            std::cout << "  Output  V(" << result.output_node
                      << (result.reference_node == "0" ? ""
                                                     : "," + result.reference_node)
                      << ") = " << result.stats.output_noise_rms << " V rms  ("
                      << result.stats.output_noise_v2 << " V^2)" << std::endl;
            if (!result.input_source.empty()) {
                std::cout << "  Input   referred to " << result.input_source
                          << "   = " << result.stats.input_noise_rms
                          << " V rms  (" << result.stats.input_noise_v2
                          << " V^2)" << std::endl;
            }
            std::cout << "  Peak of the output PSD = "
                      << std::sqrt(std::max(result.stats.peak_output_psd, 0.0))
                      << " V/sqrt(Hz) at " << result.stats.peak_output_freq
                      << " Hz" << std::endl;

            // Noise budget: the generators that actually move the number.
            if (!result.sources.empty() &&
                !result.stats.source_noise_v2.empty()) {
                std::vector<size_t> order(result.sources.size());
                for (size_t i = 0; i < order.size(); ++i) order[i] = i;
                std::sort(order.begin(), order.end(),
                          [&](size_t a, size_t b) {
                              return result.stats.source_noise_v2[a] >
                                     result.stats.source_noise_v2[b];
                          });
                const double total = std::max(result.stats.output_noise_v2, 0.0);
                std::cout << "\nNoise sources (integrated output contribution):"
                          << std::endl;
                const size_t shown = std::min<size_t>(10, order.size());
                for (size_t i = 0; i < shown; ++i) {
                    const size_t k = order[i];
                    const double v2 = result.stats.source_noise_v2[k];
                    const double pct = (total > 0.0) ? 100.0 * v2 / total : 0.0;
                    std::cout << "  " << result.sources[k].device << " ("
                              << result.sources[k].label << "): "
                              << std::sqrt(std::max(v2, 0.0)) << " V rms, "
                              << pct << "%" << std::endl;
                }
                if (order.size() > shown) {
                    std::cout << "  (" << order.size() - shown
                              << " more source(s) in the CSV)" << std::endl;
                }
            }

            const size_t nstride = std::max<size_t>(
                1, (result.points.size() + 9) / 10);
            std::cout << "\nNoise spectral density (V/sqrt(Hz)):" << std::endl;
            for (size_t i = 0; i < result.points.size(); i += nstride) {
                const auto& p = result.points[i];
                std::cout << "  f=" << p.frequency
                          << " Hz: onoise=" << p.onoise_amp
                          << " inoise=" << p.inoise_amp
                          << " |gain|=" << p.gain_mag << std::endl;
            }
            std::cout << "\nThe CSV carries the power spectral densities in "
                         "V^2/Hz (plus per-source columns when the .NOISE "
                         "line has pts_per_summary)." << std::endl;
        }
    }

    if (do_tf) {
        TFAnalysis tf(circuit);
        tf.setOutputVar(params.tf_out_var);
        tf.setInputSource(params.tf_input_source);
        tf.setVerbose(mpi.isMaster());
        if (params.opt_max_iterations > 0) tf.setMaxIterations(params.opt_max_iterations);
        TFResult result = tf.run();
        if (mpi.isMaster()) {
            if (!result.success) {
                std::cerr << "TF analysis failed: " << result.message << std::endl;
            } else {
                std::cout << "\nTransfer function analysis completed: "
                          << result.message << std::endl;
                std::cout << "  transfer function = " << result.transfer_function
                          << std::endl;
                std::cout << "  input  resistance = " << result.input_resistance
                          << " ohm" << std::endl;
                std::cout << "  output resistance = " << result.output_resistance
                          << " ohm" << std::endl;
            }
        }
    }

    mpi.barrier();
    std::cout << "[DEBUG] Rank " << mpi.getRank() << " chegou ao fim de simulateNetlist. Inciando destrutores locais...\n";

    return true;
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    if (mpi.isMaster()) {
        std::cout << "\n|--------------------------------------|\n";
        std::cout << "|           -- e-XSpice --               |\n";
        std::cout << "|   MPI + OpenMP Parallel Version        |\n";
        std::cout << "|----------------------------------------|\n\n";
    }

    mpi.printInfo();
    mpi.barrier();

    bool use_mpi = (mpi.getSize() > 1);

    if (argc >= 2 && std::string(argv[1]) == "file") {
        if (argc < 3) {
            if (mpi.isMaster()) std::cerr << "Usage: " << argv[0] << " file <netlist.cir> [output.csv]" << std::endl;
            mpi.finalize();
            return 1;
        }
        const std::string csv_path = (argc >= 4) ? argv[3] : "simulation_output.csv";
        bool success = simulateNetlist(argv[2], use_mpi, csv_path);

        std::cout << "[DEBUG] Rank " << mpi.getRank() << " sobreviveu aos destrutores. Chamando MPI_Finalize...\n";

        mpi.finalize();

        std::cout << "[DEBUG] Rank " << mpi.getRank() << " encerrou MPI com sucesso.\n";

        return success ? 0 : 1;
    }

    if (mpi.isMaster()) {
        std::cerr << "Usage: " << argv[0] << " file <netlist.cir> [output.csv]" << std::endl;
    }

    mpi.finalize();

    return 1;
}
