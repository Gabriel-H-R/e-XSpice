/**
 * @author Gabriel Henrique Silva
 * @file main.cpp
 * @brief Main entry point of the circuit simulator
 * * Supports:
 * - Programmatic circuit creation
 * - SPICE netlist parsing (.cir)
 * - DC Operating Point analysis
 * - Transient analysis
 * - MPI + OpenMP parallelization
 * * Usage:
 * ./e_xspice <mode> [options]
 * * Modes:
 * 1 - Programmatic example: Voltage divider
 * 2 - Programmatic example: Multi-Loop Network
 * 3 - Programmatic example: RC Step Response
 * 4 - Programmatic example: RLC Series Circuit
 * file <netlist.cir> - Loads and simulates netlist
 */

#include "core/Circuit.h"
#include "core/Node.h"
#include "core/NetlistParser.h"
#include "devices/Resistor.h"
#include "devices/VoltageSource.h"
#include "devices/CurrentSource.h"
#include "devices/Capacitor.h"
#include "devices/Inductor.h"
#include "analysis/DCAnalysis.h"
#include "analysis/TransientAnalysis.h"
#include "parallel/MPIManager.h"

#include <stdio.h>
#include <iostream>
#include <memory>
#include <fstream>
#include <string>

using namespace e_XSpice;

/**
 * @brief 1. Voltage Divider Circuit (divider.cir)
 */
Circuit* createDividerCircuit() {
    Circuit* circuit = new Circuit("Voltage Divider (.cir)");
    
    Node* gnd = circuit->getGroundNode();
    Node* n1 = circuit->addNode(1, "1");
    Node* n2 = circuit->addNode(2, "2");
    
    circuit->addDevice(new VoltageSource("V1", n1, gnd, 10.0));
    circuit->addDevice(new Resistor("R1", n1, n2, 1000.0));
    circuit->addDevice(new Resistor("R2", n2, gnd, 2000.0));
    
    return circuit;
}

/**
 * @brief 2. Multi-Loop Circuit (multiloop.cir)
 */
Circuit* createMultiLoopCircuit() {
    Circuit* circuit = new Circuit("Multi-Loop Circuit (.cir)");
    
    Node* gnd = circuit->getGroundNode();
    Node* n1 = circuit->addNode(1, "1");
    Node* n2 = circuit->addNode(2, "2");
    Node* n3 = circuit->addNode(3, "3");
    Node* n4 = circuit->addNode(4, "4");
    
    circuit->addDevice(new VoltageSource("V1", n1, gnd, 12.0));
    circuit->addDevice(new Resistor("R1", n1, n2, 1000.0));
    circuit->addDevice(new Resistor("R2", n2, n3, 2000.0));
    circuit->addDevice(new Resistor("R3", n3, gnd, 1500.0));
    circuit->addDevice(new Resistor("R4", n2, n4, 3000.0));
    circuit->addDevice(new Resistor("R5", n4, gnd, 2500.0));
    circuit->addDevice(new Resistor("R6", n3, n4, 1000.0));
    
    circuit->addDevice(new CurrentSource("I1", gnd, n4, 1e-3)); 
    
    return circuit;
}

/**
 * @brief 3. RC Circuit - Step Response (rc_step.cir)
 */
Circuit* createRCStepCircuit() {
    Circuit* circuit = new Circuit("RC Step Response (.cir)");
    
    Node* gnd = circuit->getGroundNode();
    Node* n1 = circuit->addNode(1, "1");
    Node* n2 = circuit->addNode(2, "2");
    
    circuit->addDevice(new VoltageSource("V1", n1, gnd, 5.0));
    circuit->addDevice(new Resistor("R1", n1, n2, 1000.0));
    circuit->addDevice(new Capacitor("C1", n2, gnd, 1e-6));
    
    return circuit;
}

/**
 * @brief 4. RLC Series Circuit (rlc_series.cir)
 */
Circuit* createRLCSeriesCircuit() {
    Circuit* circuit = new Circuit("RLC Series Circuit (.cir)");
    
    Node* gnd = circuit->getGroundNode();
    Node* n1 = circuit->addNode(1, "1");
    Node* n2 = circuit->addNode(2, "2");
    Node* n3 = circuit->addNode(3, "3");
    
    circuit->addDevice(new VoltageSource("V1", n1, gnd, 10.0));
    circuit->addDevice(new Resistor("R1", n1, n2, 100.0));
    circuit->addDevice(new Inductor("L1", n2, n3, 10e-3));
    circuit->addDevice(new Capacitor("C1", n3, gnd, 1e-6));
    
    return circuit;
}

/**
 * @brief Simulates circuit from netlist file
 */
bool simulateNetlist(const std::string& filename, bool use_mpi) {
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
    
    if (params.type == AnalysisType::DC_OP) {
        DCAnalysis dc(circuit);
        dc.setVerbose(mpi.isMaster());
        dc.setTolerance(1e-10);
        
        DCResult result = dc.run(use_mpi, false);
        
        if (mpi.isMaster()) {
            if (result.success) circuit.printSolution();
            else std::cerr << "DC Analysis failed: " << result.message << std::endl;
        }
    } else if (params.type == AnalysisType::TRANSIENT) {
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
        
        if (mpi.isMaster()) {
            if (result.success) {
                std::cout << "\nTransient analysis completed!" << std::endl;
                std::cout << "Total simulation time: " << result.total_time << " s" << std::endl;
                
                int sample_interval = std::max(1, (int)result.data.size() / 10);
                std::cout << "\nSample data points:" << std::endl;
                for (size_t i = 0; i < result.data.size(); i += sample_interval) {
                    const auto& point = result.data[i];
                    std::cout << "  t=" << point.time << "s: ";
                    for (size_t j = 0; j < std::min(size_t(3), point.node_voltages.size()); ++j) {
                        std::cout << "V" << (j+1) << "=" << point.node_voltages[j] << "V ";
                    }
                    std::cout << std::endl;
                }
            } else {
                std::cerr << "Transient Analysis failed: " << result.message << std::endl;
            }
        }
    } else {
        if (mpi.isMaster()) std::cerr << "No analysis specified in netlist (use .OP or .TRAN)" << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);
    
    if (mpi.isMaster()) {
        std::cout << "\n╔════════════════════════════════════════╗\n";
        std::cout << "║           -- e-XSpice --               ║\n";
        std::cout << "║   MPI + OpenMP Parallel Version        ║\n";
        std::cout << "╚════════════════════════════════════════╝\n\n";
    }
    
    mpi.printInfo();
    mpi.barrier();
    
    bool use_mpi = (mpi.getSize() > 1);
    
    // File mode
    if (argc >= 2 && std::string(argv[1]) == "file") {
        if (argc < 3) {
            if (mpi.isMaster()) std::cerr << "Usage: " << argv[0] << " file <netlist.cir>" << std::endl;
            mpi.finalize();
            return 1;
        }
        bool success = simulateNetlist(argv[2], use_mpi);
        mpi.finalize();
        return success ? 0 : 1;
    }
    
    // Programmatic mode
    int circuit_choice = 1;
    if (argc > 1) circuit_choice = std::atoi(argv[1]);
    
    Circuit* circuit = nullptr;
    
    switch (circuit_choice) {
        case 1:
            if (mpi.isMaster()) std::cout << "\nRunning Example 1: Voltage Divider\n" << std::endl;
            circuit = createDividerCircuit();
            break;
        case 2:
            if (mpi.isMaster()) std::cout << "\nRunning Example 2: Multi-Loop Network\n" << std::endl;
            circuit = createMultiLoopCircuit();
            break;
        case 3:
            if (mpi.isMaster()) std::cout << "\nRunning Example 3: RC Step Response\n" << std::endl;
            circuit = createRCStepCircuit();
            break;
        case 4:
            if (mpi.isMaster()) std::cout << "\nRunning Example 4: RLC Series Circuit\n" << std::endl;
            circuit = createRLCSeriesCircuit();
            break;
        default:
            if (mpi.isMaster()) std::cout << "\nInvalid choice. Defaulting to Example 1.\n" << std::endl;
            circuit_choice = 1;
            circuit = createDividerCircuit();
            break;
    }
    
    if (mpi.isMaster()) circuit->printTopology();
    
    bool sim_success = false;

    // Route logic based on analysis type needed
    if (circuit_choice == 1 || circuit_choice == 2) {
        // --- DC ANALYSIS (Static Circuits) ---
        DCAnalysis dc(*circuit);
        dc.setSolverType(SolverType::BICGSTAB);
        dc.setMaxIterations(10000);
        dc.setTolerance(1e-10);
        dc.setVerbose(mpi.isMaster());
        
        DCResult result = dc.run(use_mpi, false);
        sim_success = result.success;
        
        if (mpi.isMaster()) {
            if (sim_success) {
                circuit->printSolution();
                std::cout << "\nPerformance Summary:\n";
                std::cout << "  Assembly time: " << result.assembly_time << " s\n";
                std::cout << "  Solve time: " << result.solve_time << " s\n";
                std::cout << "  Total time: " << (result.assembly_time + result.solve_time) << " s\n";
                std::cout << "  Solver iterations: " << result.solver_result.iterations << "\n";
                std::cout << "  Final residual: " << result.solver_result.residual_norm << "\n";
            } else {
                std::cerr << "ERROR: " << result.message << std::endl;
            }
        }
    } else {
        // --- TRANSIENT ANALYSIS (Time-dependent Circuits) ---
        TransientAnalysis tran(*circuit);
        tran.setVerbose(mpi.isMaster());
        tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
        tran.setTolerance(1e-9);
        tran.setSaveInterval(10);
        
        if (mpi.isMaster()) {
            tran.setProgressCallback([](double t, int step, int total) {
                if (step % 100 == 0 || step == total) {
                    std::cout << "  Progress: " << (100.0 * step / total) << "% (t=" << t << "s)" << std::endl;
                }
            });
        }

        // Set times based on the circuit (Matching the .cir files)
        double tstart = 0.0;
        double tstep = 10e-6; // 10us
        double tstop = (circuit_choice == 3) ? 5e-3 : 10e-3; // 5ms for RC, 10ms for RLC

        TransientResult result = tran.run(tstart, tstop, tstep, use_mpi);
        sim_success = result.success;
        
        if (mpi.isMaster()) {
            if (sim_success) {
                std::cout << "\nTransient analysis completed!\n";
                std::cout << "Total simulation time: " << result.total_time << " s\n";
                
                int sample_interval = std::max(1, (int)result.data.size() / 10);
                std::cout << "\nSample data points:\n";
                for (size_t i = 0; i < result.data.size(); i += sample_interval) {
                    const auto& point = result.data[i];
                    std::cout << "  t=" << point.time << "s: ";
                    for (size_t j = 0; j < std::min(size_t(3), point.node_voltages.size()); ++j) {
                        std::cout << "V" << (j+1) << "=" << point.node_voltages[j] << "V ";
                    }
                    std::cout << "\n";
                }
            } else {
                std::cerr << "Transient Analysis failed: " << result.message << std::endl;
            }
        }
    }
    
    delete circuit;
    
    if (mpi.isMaster()) {
        std::cout << "\nSimulation completed!\n" << std::endl;
        std::cout << "Press ENTER to close..." << std::endl;
        std::cin.get();
    }
    
    mpi.finalize();
    
    return sim_success ? 0 : 1;
}