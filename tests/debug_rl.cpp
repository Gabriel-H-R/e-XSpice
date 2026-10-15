// Debug: RL step, first steps.
#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Inductor.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/analysis/TransientAnalysis.h"
#include <cstdio>

using namespace e_XSpice;

int main() {
    Circuit c("rl");
    Node* gnd = c.getGroundNode();
    Node* n1 = c.addNode(1, "1");
    Node* n2 = c.addNode(2, "2");
    c.addDevice(make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(10.0)));
    c.addDevice(make_rcp<Resistor>("R1", n1, n2, 100.0));
    c.addDevice(make_rcp<Inductor>("L1", n2, gnd, 10e-3, 0.0));

    TransientAnalysis tran(c);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    tran.setVerbose(true);

    TransientResult r = tran.run(0.0, 2.5e-5, 1e-8);
    printf("success=%d msg=%s\n", (int)r.success, r.message.c_str());
    printf("%zu points\n", r.data.size());
    for (size_t i = 0; i < r.data.size() && i < 12; ++i) {
        printf("t=%.4g  v2=%.8g  iL=%.8g  iV=%.8g\n",
               r.data[i].time, r.data[i].node_voltages[1],
               r.data[i].branch_currents[1], r.data[i].branch_currents[0]);
    }
    printf("...\n");
    for (size_t i = r.data.size() / 2; i < r.data.size(); i += r.data.size() / 8) {
        printf("t=%.4g  v2=%.8g  iL=%.8g\n",
               r.data[i].time, r.data[i].node_voltages[1], r.data[i].branch_currents[1]);
    }
    printf("final: t=%.4g iL=%.8g (expect 0.1)\n",
           r.data.back().time, r.data.back().branch_currents[1]);
    return 0;
}
