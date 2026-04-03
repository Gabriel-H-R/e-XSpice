/**
 * @file test_divider.cpp
 * @brief Unit test for the Voltage Divider circuit (.OP analysis)
 */

#include "../src/core/Circuit.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/parallel/MPIManager.h"
#include <iostream>

using namespace e_XSpice;

int main(int argc, char** argv) {
    // Initialize MPI environment
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    // Build the Voltage Divider circuit
    Circuit circuit("Voltage Divider Test");
    Node* gnd = circuit.getGroundNode();
    Node* n1 = circuit.addNode(1, "1");
    Node* n2 = circuit.addNode(2, "2");

    // V1 1 0 DC 10
    circuit.addDevice(new VoltageSource("V1", n1, gnd, 10.0));
    // R1 1 2 1k
    circuit.addDevice(new Resistor("R1", n1, n2, 1000.0));
    // R2 2 0 2k
    circuit.addDevice(new Resistor("R2", n2, gnd, 2000.0));

    // Run DC Analysis
    DCAnalysis dc(circuit);
    bool use_mpi = (mpi.getSize() > 1);
    DCResult result = dc.run(use_mpi, false);

    mpi.finalize();

    if (result.success) {
        std::cout << "Voltage Divider test PASSED!" << std::endl;
        return 0;
    } else {
        std::cerr << "Voltage Divider test FAILED: " << result.message << std::endl;
        return 1;
    }
}