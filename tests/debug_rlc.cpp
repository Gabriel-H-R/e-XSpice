// Debug: RLC response shape + pulse-edge failure.
#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Capacitor.h"
#include "../src/devices/Inductor.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/analysis/TransientAnalysis.h"
#include <cstdio>

using namespace e_XSpice;

int main() {
    Circuit c("rlc");
    Node* gnd = c.getGroundNode();
    Node* n1 = c.addNode(1, "1");
    Node* n2 = c.addNode(2, "2");
    Node* n3 = c.addNode(3, "3");
    c.addDevice(make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(10.0)));
    c.addDevice(make_rcp<Resistor>("R1", n1, n2, 100.0));
    c.addDevice(make_rcp<Inductor>("L1", n2, n3, 10e-3, 0.0, false));
    c.addDevice(make_rcp<Capacitor>("C1", n3, gnd, 1e-6, 0.0, true));

    TransientAnalysis tran(c);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    tran.setVerbose(false);
    TransientResult r = tran.run(0.0, 2e-3, 2e-7);
    printf("success=%d msg=%s\n", (int)r.success, r.message.c_str());
    printf("t(us)   v3      iL\n");
    for (size_t i = 0; i < r.data.size(); i += 200) {
        printf("%7.1f %8.4f %8.5f\n", r.data[i].time * 1e6,
               r.data[i].node_voltages[2], r.data[i].branch_currents[1]);
    }
    return 0;
}
