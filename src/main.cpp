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
#include "analysis/TransientAnalysis.h"
#include "parallel/MPIManager.h"

#include <stdio.h>
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
    } else if (!do_op && !do_dc && !do_ac) {
        if (mpi.isMaster()) std::cerr << "No analysis specified in netlist (use .OP, .TRAN, .DC or .AC)" << std::endl;
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

            const int stride = std::max<int>(1, (int)result.points.size() / 10);
            std::cout << "\nMagnitude / phase:" << std::endl;
            for (size_t i = 0; i < result.points.size(); i += stride) {
                const auto& p = result.points[i];
                std::cout << "  f=" << p.frequency << " Hz: ";
                for (size_t j = 0; j < std::min<size_t>(3, p.mag.size()); ++j) {
                    std::cout << "|" << labels[j] << "|=" << p.mag[j] << "V";
                    std::cout << "@" << p.phase_deg[j] << "deg ";
                }
                if (!p.branch_mag.empty()) {
                    std::cout << "|I|=" << p.branch_mag[0]
                              << "A@" << p.branch_phase_deg[0] << "deg ";
                }
                std::cout << std::endl;
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
