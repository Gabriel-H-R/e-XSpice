/**
 * @file test_multiloop.cpp
 * @brief Unit test for the Multi-Loop circuit (.OP analysis)
 */

#include "../src/core/Circuit.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/CurrentSource.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/parallel/MPIManager.h"
#include <iostream>

using namespace e_XSpice;

int main(int argc, char** argv) {
    // Initialize MPI environment
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    // Build the Multi-Loop circuit
    Circuit circuit("Multi-Loop Test");
    Node* gnd = circuit.getGroundNode();
    Node* n1 = circuit.addNode(1, "1");
    Node* n2 = circuit.addNode(2, "2");
    Node* n3 = circuit.addNode(3, "3");
    Node* n4 = circuit.addNode(4, "4");

    circuit.addDevice(new VoltageSource("V1", n1, gnd, 12.0));
    circuit.addDevice(new Resistor("R1", n1, n2, 1000.0));
    circuit.addDevice(new Resistor("R2", n2, n3, 2000.0));
    circuit.addDevice(new Resistor("R3", n3, gnd, 1500.0));
    circuit.addDevice(new Resistor("R4", n2, n4, 3000.0));
    circuit.addDevice(new Resistor("R5", n4, gnd, 2500.0));
    circuit.addDevice(new Resistor("R6", n3, n4, 1000.0));
    
    // I1 0 4 DC 1m
    circuit.addDevice(new CurrentSource("I1", gnd, n4, 1e-3));

    // Run DC Analysis
    DCAnalysis dc(circuit);
    bool use_mpi = (mpi.getSize() > 1);
    DCResult result = dc.run(use_mpi, false);

    mpi.finalize();

    if (result.success) {
        std::cout << "Multi-Loop test PASSED!" << std::endl;
        return 0;
    } else {
        std::cerr << "Multi-Loop test FAILED: " << result.message << std::endl;
        return 1;
    }
}