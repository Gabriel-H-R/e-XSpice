/**
 * @author Gabriel Henrique Silva
 * @file main.cpp
 * @brief Main entry point of the circuit simulator
 * 
 * Supports:
 * - Programmatic circuit creation
 * - SPICE netlist parsing (.cir)
 * - DC Operating Point analysis
 * - Transient analysis
 * - MPI + OpenMP parallelization
 * 
 * Usage:
 *   ./e_xspice <mode> [options]
 *   
 * Modes:
 *   1 - Programmatic example: Voltage divider
 *   2 - Programmatic example: Current source
 *   3 - Programmatic example: Large circuit
 *   file <netlist.cir> - Loads and simulates netlist
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

using namespace e_XSpice;

/**
 * @brief Creates example circuit 1: Voltage divider
 */
Circuit* createDividerCircuit() {
    Circuit* circuit = new Circuit("Voltage Divider");
    
    // Nodes
    Node* gnd = circuit->getGroundNode();
    Node* n1 = circuit->addNode(1, "N1");
    Node* vdd = circuit->addNode(2, "VDD");
    
    // Devices
    VoltageSource* v1 = new VoltageSource("V1", vdd, gnd, 10.0);  // 10V
    Resistor* r1 = new Resistor("R1", vdd, n1, 1000.0);           // 1kΩ
    Resistor* r2 = new Resistor("R2", n1, gnd, 2000.0);           // 2kΩ
    
    circuit->addDevice(v1);
    circuit->addDevice(r1);
    circuit->addDevice(r2);
    
    return circuit;
}

/**
 * @brief Creates example circuit 2: Circuit with current source
 */
Circuit* createCurrentSourceCircuit() {
    Circuit* circuit = new Circuit("Current Source Test");
    
    Node* gnd = circuit->getGroundNode();
    Node* n1 = circuit->addNode(1, "OUT");
    
    // 1mA flowing to the node through 10kΩ resistor
    // V = I × R = 0.001 × 10000 = 10V
    CurrentSource* i1 = new CurrentSource("I1", n1, gnd, 0.001);  // 1mA
    Resistor* r1 = new Resistor("R1", n1, gnd, 10000.0);          // 10kΩ
    
    circuit->addDevice(i1);
    circuit->addDevice(r1);
    
    return circuit;
}

/**
 * @brief Creates larger circuit for scalability testing
 */
Circuit* createLargeCircuit(int size) {
    Circuit* circuit = new Circuit("Large Resistor Network");
    
    Node* gnd = circuit->getGroundNode();
    
    // Voltage source
    Node* vdd = circuit->addNode(1, "VDD");
    VoltageSource* v1 = new VoltageSource("V1", vdd, gnd, 5.0);
    circuit->addDevice(v1);
    
    // Chain of resistors
    Node* prev = vdd;
    for (int i = 0; i < size; ++i) {
        Node* next = circuit->addNode(i + 2, "N" + std::to_string(i + 2));
        Resistor* r = new Resistor("R" + std::to_string(i), prev, next, 1000.0);
        circuit->addDevice(r);
        prev = next;
    }
    
    // Last resistor to ground
    Resistor* r_last = new Resistor("R_LAST", prev, gnd, 1000.0);
    circuit->addDevice(r_last);
    
    return circuit;
}

/**
 * @brief Simulates circuit from netlist file
 */
bool simulateNetlist(const std::string& filename, bool use_mpi) {
    MPIManager& mpi = MPIManager::getInstance();
    
    // Parse netlist
    Circuit circuit("Netlist Circuit");
    NetlistParser parser;
    parser.setVerbose(mpi.isMaster());
    
    if (mpi.isMaster()) {
        std::cout << "\nParsing netlist: " << filename << std::endl;
    }
    
    if (!parser.parse(filename, circuit)) {
        if (mpi.isMaster()) {
            std::cerr << "Parse error: " << parser.getErrorMessage() << std::endl;
        }
        return false;
    }
    
    // Prints topology
    if (mpi.isMaster()) {
        circuit.printTopology();
    }
    
    // Gets analysis type
    const auto& params = parser.getAnalysisParameters();
    
    if (params.type == AnalysisType::DC_OP) {
        // DC Analysis
        DCAnalysis dc(circuit);
        dc.setVerbose(mpi.isMaster());
        dc.setTolerance(1e-10);
        
        DCResult result = dc.run(use_mpi, false);
        
        if (mpi.isMaster()) {
            if (result.success) {
                circuit.printSolution();
            } else {
                std::cerr << "DC Analysis failed: " << result.message << std::endl;
                return false;
            }
        }
        
    } else if (params.type == AnalysisType::TRANSIENT) {
        // Transient Analysis
        TransientAnalysis tran(circuit);
        tran.setVerbose(mpi.isMaster());
        tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
        tran.setTolerance(1e-9);
        tran.setSaveInterval(10);  // Saves every 10 steps
        
        // Progress callback
        if (mpi.isMaster()) {
            tran.setProgressCallback([](double t, int step, int total) {
                if (step % 1000 == 0) {
                    double progress = 100.0 * step / total;
                    std::cout << "  Progress: " << progress << "% (t=" << t << "s)" << std::endl;
                }
            });
        }
        
        TransientResult result = tran.run(
            params.tran_tstart,
            params.tran_tstop,
            params.tran_tstep,
            use_mpi
        );
        
        if (mpi.isMaster()) {
            if (result.success) {
                std::cout << "\nTransient analysis completed!" << std::endl;
                std::cout << "Total simulation time: " << result.total_time << " s" << std::endl;
                std::cout << "Data points saved: " << result.data.size() << std::endl;
                
                // Prints some data points
                std::cout << "\nSample data points:" << std::endl;
                int sample_interval = std::max(1, (int)result.data.size() / 10);
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
                return false;
            }
        }
        
    } else {
        if (mpi.isMaster()) {
            std::cerr << "No analysis specified in netlist (use .OP or .TRAN)" << std::endl;
        }
        return false;
    }
    
    return true;
}

int main(int argc, char** argv) {
    // Initializes MPI
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);
    
    if (mpi.isMaster()) {
        std::cout << "\n";
        std::cout << "╔════════════════════════════════════════╗\n";
        std::cout << "║           -- e-XSpice --               ║\n";
        std::cout << "║   MPI + OpenMP Parallel Version        ║\n";
        std::cout << "╚════════════════════════════════════════╝\n";
        std::cout << std::endl;
    }
    
    mpi.printInfo();
    mpi.barrier();
    
    // Determines operation mode
    bool use_mpi = (mpi.getSize() > 1);
    
    if (argc >= 2 && std::string(argv[1]) == "file") {
        // Netlist mode
        if (argc < 3) {
            if (mpi.isMaster()) {
                std::cerr << "Usage: " << argv[0] << " file <netlist.cir>" << std::endl;
            }
            mpi.finalize();
            return 1;
        }
        
        std::string filename = argv[2];
        bool success = simulateNetlist(filename, use_mpi);
        
        mpi.finalize();
        return success ? 0 : 1;
    }
    
    // Programmatic examples mode (original)
    int circuit_choice = 1;
    if (argc > 1) {
        circuit_choice = std::atoi(argv[1]);
    }
    
    Circuit* circuit = nullptr;
    
    if (circuit_choice == 1) {
        if (mpi.isMaster()) {
            std::cout << "\nRunning Example 1: Voltage Divider\n" << std::endl;
        }
        circuit = createDividerCircuit();
    } else if (circuit_choice == 2) {
        if (mpi.isMaster()) {
            std::cout << "\nRunning Example 2: Current Source\n" << std::endl;
        }
        circuit = createCurrentSourceCircuit();
    } else {
        int size = 100;
        if (argc > 2) {
            size = std::atoi(argv[2]);
        }
        if (mpi.isMaster()) {
            std::cout << "\nRunning Example 3: Large Circuit (" << size << " resistors)\n" << std::endl;
        }
        circuit = createLargeCircuit(size);
    }
    
    // Prints topology (master only)
    if (mpi.isMaster()) {
        circuit->printTopology();
    }
    
    // Creates DC analysis
    DCAnalysis dc(*circuit);
    dc.setSolverType(SolverType::BICGSTAB);
    dc.setMaxIterations(10000);
    dc.setTolerance(1e-10);
    dc.setVerbose(mpi.isMaster());  // Verbose only on master
    
    // Runs analysis
    use_mpi = (mpi.getSize() > 1);
    DCResult result = dc.run(use_mpi, false);
    
    // Prints results
    if (mpi.isMaster()) {
        if (result.success) {
            circuit->printSolution();
            
            std::cout << "\nPerformance Summary:" << std::endl;
            std::cout << "  Assembly time: " << result.assembly_time << " s" << std::endl;
            std::cout << "  Solve time: " << result.solve_time << " s" << std::endl;
            std::cout << "  Total time: " << (result.assembly_time + result.solve_time) << " s" << std::endl;
            std::cout << "  Solver iterations: " << result.solver_result.iterations << std::endl;
            std::cout << "  Final residual: " << result.solver_result.residual_norm << std::endl;
        } else {
            std::cerr << "ERROR: " << result.message << std::endl;
        }
    }
    
    // Cleanup
    delete circuit;
    
    if (mpi.isMaster()) {
        std::cout << "\nSimulation completed successfully!\n" << std::endl;

        std::cout << "Press ENTER to close..." << std::endl;
        std::cin.get();
    }
    
    mpi.finalize();
    
    return result.success ? 0 : 1;
}
