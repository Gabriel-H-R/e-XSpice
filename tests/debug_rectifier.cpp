// Debug: the netlist rectifier with the pulse edge that fails to converge.
#include "../src/core/Circuit.h"
#include "../src/core/NetlistParser.h"
#include "../src/devices/Diode.h"
#include "../src/analysis/TransientAnalysis.h"
#include <cstdio>

using namespace e_XSpice;

int main() {
    Diode::setDebug(false);

    const std::string cir = R"(
.MODEL DMOD D (IS=1e-12 N=1.0 CJO=1e-9 VJ=0.7)
V1 1 0 PULSE(0 5 1u 10u 10u 500u 1m)
R1 1 2 1k
D1 2 3 DMOD
C1 3 0 100n IC=0
L1 3 4 1m IC=0
RL 4 0 100
.TRAN 1u 2m 0
.END
)";
    Circuit c("netlist");
    NetlistParser p;
    if (!p.parseString(cir, c)) {
        printf("parse error: %s\n", p.getErrorMessage().c_str());
        return 1;
    }

    TransientAnalysis tran(c);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    tran.setVerbose(true);
    TransientResult r = tran.run(0.0, 2e-3, 1e-6);
    printf("success=%d msg=%s\n", (int)r.success, r.message.c_str());
    if (!r.data.empty()) {
        printf("t(us)  v2      v3      iV\n");
        for (size_t i = 0; i < r.data.size(); i += 50) {
            printf("%7.1f %7.4f %7.4f %9.6f\n", r.data[i].time * 1e6,
                   r.data[i].node_voltages[1], r.data[i].node_voltages[2],
                   r.data[i].branch_currents[0]);
        }
    }
    return 0;
}
