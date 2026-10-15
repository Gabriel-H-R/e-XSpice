// Debug: DC convergence of a current-source-biased BJT.
#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/CurrentSource.h"
#include "../src/devices/Bjt.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/models/core_models/BjtModel.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/parallel/MPIManager.h"
#include <cstdio>

using namespace e_XSpice;

int main(int argc, char** argv) {
    MPIManager::getInstance().initialize(&argc, &argv);
    Bjt::setDebug(true);

    auto model = make_rcp<BjtModel>("QN", false);
    std::string msg;
    if (!model->processParams(msg)) {
        printf("model error: %s\n", msg.c_str());
        return 1;
    }
    model->BF = 100.0;
    model->IS = 1e-16;
    model->processParams(msg);

    Circuit c("ce");
    Node* gnd = c.getGroundNode();
    Node* nrail = c.addNode(1, "1");
    Node* ncol = c.addNode(2, "2");
    Node* nbase = c.addNode(3, "3");
    c.addDevice(make_rcp<VoltageSource>("V1", nrail, gnd, SourceWaveform::dc(10.0)));
    c.addDevice(make_rcp<Resistor>("RC", nrail, ncol, 1000.0));
    c.addDevice(make_rcp<CurrentSource>("IB", gnd, nbase, SourceWaveform::dc(10e-6)));
    c.addDevice(make_rcp<Bjt>("Q1", ncol, nbase, gnd, model));

    DCAnalysis dc(c);
    dc.setVerbose(true);
    DCResult r = dc.run(false, false);
    printf("success=%d msg=%s\n", (int)r.success, r.message.c_str());
    printf("v2=%.9g v3=%.9g\n", ncol->getVoltage(), nbase->getVoltage());
    Bjt* q = static_cast<Bjt*>(c.getDevice("Q1"));
    printf("ib=%.9g ic=%.9g\n", q->getBaseCurrent(), q->getCollectorCurrent());
    MPIManager::getInstance().finalize();
    return 0;
}
