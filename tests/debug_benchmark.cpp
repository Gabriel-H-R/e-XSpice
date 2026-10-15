// Regression test: diode with grounded anode in transient analysis.
//
// Covers the heap-corruption bug where Diode::loadQ wrote qd to Q[-1] when
// the anode is ground (Pri aliases the anode row, which does not exist).
// The netlist below mirrors examples/benchmark_rlcd.cir (grounded-anode
// clipper D2, floating diode D1, L/C stages, SIN source).
#include "../src/core/Circuit.h"
#include "../src/core/NetlistParser.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"
#include <cmath>
#include <cstdio>

using namespace e_XSpice;

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    // Same topology as benchmark_rlcd.cir, shortened stop time.
    const std::string cir = R"(
.MODEL DMODEL D ( IS=1e-14 N=1.0 CJO=10p VJ=0.7 )
V1 1 0 SIN 0 20 1000
R1 1 2 100
L1 2 3 15m
C1 3 0 10u
R2 3 4 220
L2 4 5 5m
C2 5 0 4.7u
D1 5 6 DMODEL
D2 0 5 DMODEL
Vbias 6 0 DC 5.0
R3 5 7 1k
C3 7 0 1u
L3 7 8 1m
Rload 8 0 4.7k
.TRAN 5u 5m
.END
)";
    Circuit circuit("GroundedAnode");
    NetlistParser parser;
    if (!parser.parseString(cir, circuit)) {
        printf("FAIL: parse error: %s\n", parser.getErrorMessage().c_str());
        mpi.finalize();
        return 1;
    }

    TransientAnalysis tran(circuit);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    tran.setTolerance(1e-9);
    tran.setSaveInterval(10);
    tran.setVerbose(false);

    TransientResult result = tran.run(0.0, 5e-3, 5e-6, false);
    if (!result.success) {
        printf("FAIL: transient: %s\n", result.message.c_str());
        mpi.finalize();
        return 1;
    }
    if (result.data.empty()) {
        printf("FAIL: no data points recorded\n");
        mpi.finalize();
        return 1;
    }
    for (const auto& point : result.data) {
        for (double v : point.node_voltages) {
            if (!std::isfinite(v) || std::abs(v) > 1e6) {
                printf("FAIL: non-physical voltage %g at t=%g\n", v, point.time);
                mpi.finalize();
                return 1;
            }
        }
    }

    printf("Grounded-anode diode transient PASSED (%llu points)\n",
           (unsigned long long)result.data.size());
    mpi.finalize();
    return 0;
}
