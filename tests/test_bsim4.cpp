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
 *  5. Common-source DC transfer (NMOS and PMOS) with KCL closure: the drain
 *     current reported by the device must match what the external resistors
 *     carry. This is the check that a mis-stamped residual cannot survive
 *     (the current used to be evaluated correctly but written into the wrong
 *     rows, leaving the transistor effectively open with every "Vout is high"
 *     assertion still passing).
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
#include "../src/analysis/DCSweep.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"

#include <algorithm>
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

    // ---- 5) common-source DC transfer + KCL closure (NMOS/PMOS) ----
    //
    // Regression guard for the STAMP, not just convergence: the level-2
    // device used to report a correct drain current (609 uA at Vgs=2 V)
    // while stamping it into the wrong residual rows, so every "Vout is
    // high at Vin=0" check above still passed with an effectively open
    // transistor. Here the drain current is compared with the current the
    // external resistors must carry (KCL closure), which cannot hold unless
    // the terminal currents reached the right rows with the right sign.
    if (master) std::printf("\n5) Common-source transfer + KCL closure\n");
    {
        // VDD --Rload-- drain, gate swept; the drain is a REAL node (not
        // pinned by a source), so the operating point is set by the device.
        //   NMOS: source/bulk at ground, drain loaded from the rail.
        //   PMOS: source/bulk at the rail, gate sweep pulls the drain down.
        auto build = [](Circuit& c, bool nmos, double vg) {
            ModelRegistry reg;
            RCP<const MosfetBSIM4Model> m =
                bsimCard(reg, nmos ? "NB5" : "PB5", nmos);
            Node* gnd = c.getGroundNode();
            Node* rail = c.addNode(1, "RAIL");
            Node* g = c.addNode(2, "G");
            Node* d = c.addNode(3, "D");
            c.addDevice(make_rcp<VoltageSource>("V1", rail, gnd,
                                                 SourceWaveform::dc(3.0)));
            c.addDevice(make_rcp<VoltageSource>("V2", g, gnd,
                                                 SourceWaveform::dc(vg)));
            // NMOS: pull-up load from the rail (sink device).
            // PMOS: pull-up device with the drain loaded to ground — the
            // only wiring where the rail actually drives current through it.
            if (nmos) {
                c.addDevice(make_rcp<Resistor>("R1", rail, d, 1e3));
            } else {
                c.addDevice(make_rcp<Resistor>("R1", d, gnd, 1e3));
            }
            c.addDevice(make_rcp<Resistor>("RG", d, gnd, 1e6));
            MosfetInstanceParams in;
            in.W = 2e-6; in.L = 0.5e-6;
            c.addDevice(make_rcp<Mosfet>("M1", d, g, nmos ? gnd : rail,
                                         nmos ? gnd : rail,
                                         RCP<const DeviceModel>(m), in));
        };

        // Current the external resistors force into the drain node at Vd.
        auto expectedId = [](bool nmos, double vd) {
            return nmos ? (3.0 - vd) / 1e3 - vd / 1e6
                        : -(vd / 1e3 + vd / 1e6);
        };

        // ---- transfer curve (DCSweep): monotone Vd in the right direction
        auto transfer = [&](bool nmos, std::vector<double>& vd) {
            Circuit c("cs");
            build(c, nmos, 0.0);
            const int id = nodeIndex(c, c.getNode(3));
            DCSweepSpec spec;
            spec.source = "V2";
            spec.start = 0.0; spec.stop = 2.0; spec.step = 0.5;
            DCSweep sweep(c);
            sweep.setVerbose(false);
            sweep.setSweep(spec);
            DCSweepResult r = sweep.run(use_mpi);
            vd.clear();
            if (!r.success) return false;
            for (const auto& p : r.points) vd.push_back(p.node_voltages[id]);
            return true;
        };

        {   // NMOS: pulling the gate up must pull the drain DOWN
            std::vector<double> vd;
            const bool ran = transfer(true, vd);
            checkTrue("NMOS sweep converges", ran);
            if (ran && vd.size() >= 5) {
                bool monotone = true;
                for (size_t i = 1; i < vd.size(); ++i) {
                    if (vd[i] >= vd[i - 1]) monotone = false;
                }
                checkTrue("Vd falls monotonically with Vgs", monotone);
                checkTrue("NMOS conducts (Vd drops below 2.8 V)",
                          vd.back() < 2.8 && vd.back() < 0.9 * vd.front());
            }
            ok = ok && (g_failures == 0);
        }
        {   // PMOS pull-up with the drain loaded to ground: raising the gate
            // shrinks VSG, the device weakens and the drain falls
            std::vector<double> vd;
            const bool ran = transfer(false, vd);
            checkTrue("PMOS sweep converges", ran);
            if (ran && vd.size() >= 5) {
                bool monotone = true;
                for (size_t i = 1; i < vd.size(); ++i) {
                    if (vd[i] >= vd[i - 1]) monotone = false;
                }
                checkTrue("Vd falls monotonically with Vg (PMOS weakens)",
                          monotone);
                checkTrue("PMOS pulls the drain up at Vsg=3 V",
                          vd.front() > 2.0 * vd.back());
            }
            ok = ok && (g_failures == 0);
        }

        // ---- KCL closure at one converged operating point ----------------
        // KCL at the drain node: (3 - Vd)/Rload = Id + Vd/Rprobe. This can
        // only hold if the terminal current reached the right residual row
        // with the right polarity sign.
        auto kclCheck = [&](bool nmos, double vg, double& id_out,
                            double& vd_out) {
            Circuit c("cs");
            build(c, nmos, vg);
            DCAnalysis dc(c);
            dc.setVerbose(false);
            DCResult r = dc.run(use_mpi, false);
            if (!r.success) return false;
            const int id = nodeIndex(c, c.getNode(3));
            vd_out = r.node_voltages[id];
            auto* mos = dynamic_cast<Mosfet*>(c.getDevice("M1"));
            if (!mos) return false;
            id_out = mos->getDrainCurrent();
            return true;
        };

        {
            double idv = 0.0, vdv = 0.0;
            const bool ran = kclCheck(true, 2.0, idv, vdv);
            checkTrue("NMOS DC converges at Vgs=2 V", ran);
            if (ran) {
                const double want = expectedId(true, vdv);
                checkTrue("NMOS KCL closes (stamp + polarity)",
                          std::abs(want - idv) < 1e-9);
                if (std::abs(want - idv) >= 1e-9) {
                    std::printf("       |KCL - Id| = %.6g A (Id=%.6g)\n",
                                std::abs(want - idv), idv);
                }
                checkRange("NMOS current at Vgs=2 V", idv, 2e-4, 2e-3);
            }
            ok = ok && (g_failures == 0);
        }
        {
            double idv = 0.0, vdv = 0.0;
            const bool ran = kclCheck(false, 0.0, idv, vdv);
            checkTrue("PMOS DC converges at Vsg=3 V", ran);
            if (ran) {
                const double want = expectedId(false, vdv);
                checkTrue("PMOS KCL closes (stamp + polarity)",
                          std::abs(want - idv) < 1e-9);
                if (std::abs(want - idv) >= 1e-9) {
                    std::printf("       |KCL - Id| = %.6g A (Id=%.6g)\n",
                                std::abs(want - idv), idv);
                }
                // A PMOS sources the drain node: negative drain current.
                checkTrue("PMOS current flows out of the drain",
                          idv < -1e-4);
            }
            ok = ok && (g_failures == 0);
        }
    }

    if (master) {
        std::printf(ok && g_failures == 0 ? "BSIM4 level-2 PASSED\n"
                                          : "BSIM4 level-2 FAILED\n");
    }
    mpi.finalize();
    return (ok && g_failures == 0) ? 0 : 1;
}
