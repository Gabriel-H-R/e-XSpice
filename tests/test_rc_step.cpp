/**
 * @file test_rc_step.cpp
 * @brief Unit test for the RC Step Response circuit (.TRAN analysis)
 */

#include "../src/core/Circuit.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Capacitor.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"
#include <iostream>

using namespace e_XSpice;

int main(int argc, char** argv) {
    // Initialize MPI environment
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    // Build the RC Step circuit
    Circuit circuit("RC Step Response Test");
    Node* gnd = circuit.getGroundNode();
    Node* n1 = circuit.addNode(1, "1");
    Node* n2 = circuit.addNode(2, "2");

    circuit.addDevice(new VoltageSource("V1", n1, gnd, 5.0));
    circuit.addDevice(new Resistor("R1", n1, n2, 1000.0));
    circuit.addDevice(new Capacitor("C1", n2, gnd, 1e-6)); // 1uF

    // Run Transient Analysis
    TransientAnalysis tran(circuit);
    bool use_mpi = (mpi.getSize() > 1);
    
    // .TRAN 10u 5m 0 (tstep = 10us, tstop = 5ms, tstart = 0)
    double tstart = 0.0;
    double tstop = 5e-3;
    double tstep = 10e-6;
    
    TransientResult result = tran.run(tstart, tstop, tstep, use_mpi);

    mpi.finalize();

    if (result.success) {
        std::cout << "RC Step Response test PASSED!" << std::endl;
        return 0;
    } else {
        std::cerr << "RC Step Response test FAILED: " << result.message << std::endl;
        return 1;
    }
}