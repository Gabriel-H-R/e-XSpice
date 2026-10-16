/**
 * @file test_bsim4.cpp
 * @brief Verification of the BSIM4 level-2 MOSFET path (NMOS2/PMOS2).
 *
 * Every check is analytic or a physics invariant:
 *  1. Registry mechanics (NMOS2/PMOS2 factories, LEVEL=2 required).
 *  2. Type safety: BSIM4 Mosfet reports isBSIM4() and a null level1Model()
 *     (dynamic_cast, never static_cast across sibling card classes);
 *     level-1 Mosfet reports the opposite.
 *  3. CMOS inverter DC (Vin=0): converges, Vout pulled high.
 *  4. CMOS inverter transient smoke: converges, finite, output stays high.
 *
 * Runs under plain ctest (size==1, serial path) and under
 * `mpiexec -n 4` (exercises Mosfet::beginSolve/storeState on the BSIM4
 * branch on every rank — the former 0xC0000005 MPI crash site).
 */

#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Mosfet.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/models/core_models/MosfetModel.h"
#include "../src/models/core_models/MosfetBSIM4Model.h"
#include "../src/models/core/ModelRegistry.h"
#include "../src/models/core/ModelCard.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace e_XSpice;

static int g_failures = 0;

static void checkTrue(const char* what, bool cond) {
    if (!cond) ++g_failures;
    std::printf("  [%s] %s\n", cond ? "OK" : "FAIL", what);
}

static void checkRange(const char* what, double got, double lo, double hi) {
    const bool ok = std::isfinite(got) && got >= lo && got <= hi;
    if (!ok) ++g_failures;
    std::printf("  [%s] %-40s got=% .8g in [% .3g, % .3g]\n",
                ok ? "OK" : "FAIL", what, got, lo, hi);
}

// Minimal BSIM4 card through the registry (mirrors bsim4_test_netlist.cir).
static RCP<const MosfetBSIM4Model> bsimCard(ModelRegistry& reg,
                                            const std::string& name,
                                            bool nmos) {
    ModelCard card;
    card.name = name;
    card.type = nmos ? "NMOS2" : "PMOS2";
    card.params["LEVEL"] = 2.0;
    card.params["VTH0"] = nmos ? 0.4 : -0.4;
    card.params["KP"] = nmos ? 100e-6 : 40e-6;
    reg.createModel(card);
    return reg.findAs<MosfetBSIM4Model>(name);
}

static int nodeIndex(const Circuit& c, const Node* n) {
    const auto& nodes = c.getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].get() == n) return static_cast<int>(i);
    }
    return -1;
}

// CMOS inverter: VDD=5, VIN in, VOUT out, 1M probe to ground.
static void buildInverter(Circuit& c,
                          RCP<const MosfetBSIM4Model> nm,
                          RCP<const MosfetBSIM4Model> pm,
                          double vin, Node*& vout) {
    Node* gnd = c.getGroundNode();
    Node* vdd = c.addNode(1, "VDD");
    Node* vinN = c.addNode(2, "VIN");
    vout = c.addNode(3, "VOUT");
    c.addDevice(make_rcp<VoltageSource>("V1", vdd, gnd,
                                        SourceWaveform::dc(5.0)));
    c.addDevice(make_rcp<VoltageSource>("V2", vinN, gnd,
                                        SourceWaveform::dc(vin)));
    MosfetInstanceParams in, ip;
    in.W = 100e-6; in.L = 1e-6;
    ip.W = 200e-6; ip.L = 1e-6;
    c.addDevice(make_rcp<Mosfet>("M1", vout, vinN, gnd, gnd,
                                 RCP<const DeviceModel>(nm), in));
    c.addDevice(make_rcp<Mosfet>("M2", vout, vinN, vdd, vdd,
                                 RCP<const DeviceModel>(pm), ip));
    c.addDevice(make_rcp<Resistor>("R1", vout, gnd, 1e6));
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);
    const bool use_mpi = (mpi.getSize() > 1);
    const bool master = mpi.isMaster();
    bool ok = true;

    if (master) std::printf("BSIM4 level-2 checks (size=%d)\n", mpi.getSize());

    // ---- 1) registry ----
    {
        if (master) std::printf("\n1) Registry mechanics (NMOS2/PMOS2)\n");
        ModelRegistry reg;
        checkTrue("NMOS2 factory registered", reg.hasModelType("NMOS2"));
        checkTrue("PMOS2 factory registered", reg.hasModelType("PMOS2"));
        RCP<const MosfetBSIM4Model> nm = bsimCard(reg, "NB", true);
        RCP<const MosfetBSIM4Model> pm = bsimCard(reg, "PB", false);
        checkTrue("NMOS2 card created", (bool)nm);
        checkTrue("PMOS2 card created", (bool)pm);
        if (nm) checkTrue("card reports isBSIM4()", nm->isBSIM4());
        ok = ok && (g_failures == 0);

        bool badLevel = false;
        try {
            ModelCard c; c.name = "BL"; c.type = "NMOS2";
            c.params["LEVEL"] = 1.0;
            reg.createModel(c);
        } catch (const std::exception&) { badLevel = true; }
        checkTrue("LEVEL=1 rejected on NMOS2", badLevel);
        ok = ok && (g_failures == 0);
    }

    // ---- 2) type safety (dynamic_cast, never cross-sibling static_cast) ----
    {
        if (master) std::printf("\n2) level1Model() type safety\n");
        ModelRegistry reg;
        RCP<const MosfetBSIM4Model> nm = bsimCard(reg, "NB2", true);
        Circuit c("t");
        Node* gnd = c.getGroundNode();
        Node* d = c.addNode(1, "d");
        Node* g = c.addNode(2, "g");
        MosfetInstanceParams inst;
        auto dev = make_rcp<Mosfet>("M", d, g, gnd, gnd,
                                    RCP<const DeviceModel>(nm), inst);
        checkTrue("BSIM4 Mosfet reports isBSIM4()", dev->isBSIM4());
        checkTrue("BSIM4 level1Model() is null (dynamic_cast)",
                  !dev->level1Model());

        ModelRegistry reg1;
        ModelCard card;
        card.name = "N1"; card.type = "NMOS";
        card.params["VTO"] = 1.0; card.params["KP"] = 2e-5;
        reg1.createModel(card);
        RCP<const MosfetModel> l1 = reg1.findAs<MosfetModel>("N1");
        checkTrue("level-1 card created", (bool)l1);
        auto dev1 = make_rcp<Mosfet>("M1", d, g, gnd, gnd,
                                     RCP<const DeviceModel>(l1), inst);
        checkTrue("level-1 Mosfet is not BSIM4", !dev1->isBSIM4());
        checkTrue("level-1 level1Model() is non-null", (bool)dev1->level1Model());
        ok = ok && (g_failures == 0);
    }

    // ---- 3) inverter DC (Vin=0 -> Vout high) ----
    {
        if (master) std::printf("\n3) CMOS inverter DC (BSIM4)\n");
        ModelRegistry reg;
        RCP<const MosfetBSIM4Model> nm = bsimCard(reg, "NB3", true);
        RCP<const MosfetBSIM4Model> pm = bsimCard(reg, "PB3", false);
        Circuit c("inv");
        Node* vout = nullptr;
        buildInverter(c, nm, pm, 0.0, vout);

        DCAnalysis dc(c);
        dc.setVerbose(false);
        DCResult r = dc.run(use_mpi, false);
        checkTrue("DC converges", r.success);
        ok = ok && r.success;
        if (r.success) {
            checkRange("Vout pulled high", vout->getVoltage(), 4.5, 6.0);
            bool finite = true;
            for (const auto& n : c.getNodes()) {
                if (!std::isfinite(n->getVoltage())) { finite = false; break; }
            }
            checkTrue("all node voltages finite", finite);
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 4) inverter transient smoke ----
    {
        if (master) std::printf("\n4) CMOS inverter transient (BSIM4)\n");
        ModelRegistry reg;
        RCP<const MosfetBSIM4Model> nm = bsimCard(reg, "NB4", true);
        RCP<const MosfetBSIM4Model> pm = bsimCard(reg, "PB4", false);
        Circuit c("invtr");
        Node* vout = nullptr;
        buildInverter(c, nm, pm, 0.0, vout);
        const int idx = nodeIndex(c, vout);

        TransientAnalysis tran(c);
        tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
        tran.setTolerance(1e-9);
        tran.setSaveInterval(10);
        tran.setVerbose(false);
        TransientResult r = tran.run(0.0, 100e-9, 1e-9, use_mpi);
        checkTrue("transient converges", r.success && !r.data.empty());
        ok = ok && r.success && !r.data.empty();
        if (r.success && !r.data.empty() && idx >= 0) {
            bool finite = true;
            for (const auto& p : r.data) {
                for (double v : p.node_voltages) {
                    if (!std::isfinite(v)) { finite = false; break; }
                }
                if (!finite) break;
            }
            checkTrue("all transient voltages finite", finite);
            checkRange("Vout stays high at tstop",
                       r.data.back().node_voltages[idx], 4.0, 6.5);
        }
        ok = ok && (g_failures == 0);
    }

    if (master) {
        std::printf(ok && g_failures == 0 ? "BSIM4 level-2 PASSED\n"
                                          : "BSIM4 level-2 FAILED\n");
    }
    mpi.finalize();
    return (ok && g_failures == 0) ? 0 : 1;
}
