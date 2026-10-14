// Debug harness: divider DC problem, prints the system and iterations.
#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/analysis/DaeSystem.h"
#include <iostream>

using namespace e_XSpice;

int main() {
    Circuit circuit("Divider");
    Node* gnd = circuit.getGroundNode();
    Node* n1 = circuit.addNode(1, "1");
    Node* n2 = circuit.addNode(2, "2");

    circuit.addDevice(new VoltageSource("V1", n1, gnd, SourceWaveform::dc(10.0)));
    circuit.addDevice(new Resistor("R1", n1, n2, 1000.0));
    circuit.addDevice(new Resistor("R2", n2, gnd, 2000.0));

    DaeSystem sys(circuit);
    sys.setup();
    sys.setMaxIterations(50);
    sys.setRelTol(1e-3);
    sys.setAbsTol(1e-12);
    sys.setGminSchedule(1e-3, 1e-12);
    sys.setVerbose(true);
    std::cout << "system size = " << sys.size() << "\n";

    NLSolveResult r = sys.solveNonlinear(0.0, 0.0, false);
    std::cout << "converged=" << r.converged << " iters=" << r.total_iterations
              << " msg=" << r.message << "\n";
    {
        const auto& J = sys.jacobianF();
        std::cout << "nnz=" << J.nonZeros() << "\n";
        for (int row = 0; row < sys.size(); ++row) {
            const auto& v = J.values();
            const auto& cols = J.colIndices();
            const auto& rp = J.rowPointers();
            std::cout << "row " << row << ":";
            for (int k = rp[row]; k < rp[row + 1]; ++k) {
                std::cout << " (" << cols[k] << "," << v[k] << ")";
            }
            std::cout << "  rhs=" << sys.rhs()[row] << "\n";
        }
        // Direct check: A * x_prev vs b
        std::vector<double> check;
        J.multiplyVector(std::vector<double>(sys.size(), 0.0), check);
        std::cout << "A*0 = [" << check[0] << "," << check[1] << "," << check[2] << "]\n";
        LinearSolver ls(SolverType::BICGSTAB);
        ls.setMaxIterations(200);
        ls.setTolerance(1e-12);
        std::vector<double> sol;
        SolverResult lr = ls.solve(J, sys.rhs(), sol);
        std::cout << "direct solve: " << lr.converged << " it=" << lr.iterations
                  << " msg=" << lr.message << " sol=[" << sol[0] << ","
                  << sol[1] << "," << sol[2] << "]\n";
    }
    return 0;
}