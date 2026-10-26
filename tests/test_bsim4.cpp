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
#include <limits>
#include <string>
#include <utility>
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

    // ---- 6) P0 regression guards --------------------------------------
    //
    // These lock in the fixes of the "P0" audit batch.  Each one failed
    // before the fix:
    //   a) junction current scaled with M^2 (AREA was fed M and NF was fed M)
    //   b) U0 was ignored by the channel current (beta used SPICE3 KP)
    //   c) the DC threshold read GAMMA while Abulk read K1 -> no body effect
    //      by default, and the two blocks disagreed
    //   d) the junction depletion charge was always zero (bulkCap never
    //      populated, and CBS was multiplied by its own temperature factor)
    //   e) the RDSMOD/RGATEMOD internal nodes were never stamped (drain
    //      resistor behind `if (false)`, rg_/rd_ never set, resistances used
    //      as conductances) and the Jacobian slots were misaligned by 4
    if (master) std::printf("\n6) P0 regression guards\n");
    {
        using KV = std::pair<std::string, double>;
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        auto cardX = [](ModelRegistry& reg, const std::string& name,
                        bool nmos, const std::vector<KV>& extra) {
            ModelCard card;
            card.name = name;
            card.type = nmos ? "NMOS2" : "PMOS2";
            card.params["LEVEL"] = 2.0;
            card.params["VTH0"] = nmos ? 0.4 : -0.4;
            card.params["KP"] = nmos ? 100e-6 : 40e-6;
            for (const auto& kv : extra) card.params[kv.first] = kv.second;
            reg.createModel(card);
            return reg.findAs<MosfetBSIM4Model>(name);
        };

        // Junction-only operating point: gate far below Vth, drain biased
        // against the bulk by an ideal source, so only the bulk diode runs.
        auto junctionOp = [&](bool nmos, const MosfetInstanceParams& in,
                              const std::vector<KV>& extra, double vd,
                              double& id, double& qd) {
            ModelRegistry reg;
            RCP<const MosfetBSIM4Model> m = cardX(reg, "NJ", nmos, extra);
            Circuit c("j");
            Node* gnd = c.getGroundNode();
            Node* d = c.addNode(1, "D");
            Node* g = c.addNode(2, "G");
            c.addDevice(make_rcp<VoltageSource>("V1", d, gnd,
                                                 SourceWaveform::dc(vd)));
            c.addDevice(make_rcp<VoltageSource>("V2", g, gnd,
                                                 SourceWaveform::dc(0.0)));
            (void)0;
            c.addDevice(make_rcp<Mosfet>("M1", d, g, gnd, gnd,
                                         RCP<const DeviceModel>(m), in));
            DCAnalysis dc(c);
            dc.setVerbose(false);
            DCResult r = dc.run(false, false);
            id = NaN; qd = NaN;
            if (!r.success) return false;
            auto* mos = dynamic_cast<Mosfet*>(c.getDevice("M1"));
            if (!mos) return false;
            id = mos->getDrainCurrent();
            qd = mos->getChargeDrain();
            return true;
        };

        // Common-source stage used by the U0 / body-effect / series-R checks.
        auto commonSource = [&](const std::vector<KV>& extra,
                                const MosfetInstanceParams& in,
                                double vbs, double& vd, double& id) {
            ModelRegistry reg;
            RCP<const MosfetBSIM4Model> m = cardX(reg, "NC", true, extra);
            Circuit c("cs");
            Node* gnd = c.getGroundNode();
            Node* rail = c.addNode(1, "RAIL");
            Node* g = c.addNode(2, "G");
            Node* d = c.addNode(3, "D");
            Node* b = c.addNode(4, "B");
            c.addDevice(make_rcp<VoltageSource>("V1", rail, gnd,
                                                 SourceWaveform::dc(3.0)));
            c.addDevice(make_rcp<VoltageSource>("V2", g, gnd,
                                                 SourceWaveform::dc(2.0)));
            c.addDevice(make_rcp<VoltageSource>("V3", b, gnd,
                                                 SourceWaveform::dc(vbs)));
            c.addDevice(make_rcp<Resistor>("R1", rail, d, 1e3));
            c.addDevice(make_rcp<Resistor>("RG", d, gnd, 1e6));
            c.addDevice(make_rcp<Mosfet>("M1", d, g, gnd, b,
                                         RCP<const DeviceModel>(m), in));
            DCAnalysis dc(c);
            dc.setVerbose(false);
            DCResult r = dc.run(false, false);
            vd = NaN; id = NaN;
            if (!r.success) return false;
            const int idx = nodeIndex(c, c.getNode(3));
            vd = r.node_voltages[idx];
            auto* mos = dynamic_cast<Mosfet*>(c.getDevice("M1"));
            if (!mos) return false;
            id = mos->getDrainCurrent();
            return true;
        };

        // (a) M scales the junction current linearly, not quadratically.
        {
            const std::vector<KV> jpar = {{"VTH0", 1.5}, {"JSS", 1.0}};
            MosfetInstanceParams in1, in2;
            in1.W = 100e-6; in1.L = 1e-6; in1.M = 1.0; in1.NF = 1.0;
            in1.AD = 1e-6; in1.AS = 1e-6;
            in2 = in1; in2.M = 2.0;
            double i1 = 0.0, q1 = 0.0, i2 = 0.0, q2 = 0.0;
            const bool ok1 = junctionOp(true, in1, jpar, 1.0, i1, q1);
            const bool ok2 = junctionOp(true, in2, jpar, 1.0, i2, q2);
            checkTrue("junction-only DC converges (M=1 and M=2)", ok1 && ok2);
            checkTrue("junction current is non-zero (JSS is honoured)",
                      ok1 && std::abs(i1) > 1e-9);
            checkTrue("M=2 junction current is exactly 2x M=1 (was M^2)",
                      ok1 && ok2 &&
                      std::abs(i2 - 2.0 * i1) < 1e-3 * std::abs(i1));
            checkTrue("M=2 junction charge is exactly 2x M=1",
                      ok1 && ok2 &&
                      std::abs(q2 - 2.0 * q1) < 1e-3 * std::abs(q1));
        }

        // (b) A U0-only card must match the equivalent KP card.
        {
            const double cox = 3.9 * 8.854214871e-12 / 30.0e-10;
            const double u0 = 0.05;
            MosfetInstanceParams in;
            in.W = 2e-6; in.L = 0.5e-6;
            double vdu = 0.0, idu = 0.0, vdk = 0.0, idk = 0.0;
            const bool oku = commonSource({{"U0", u0}}, in, 0.0, vdu, idu);
            const bool okk = commonSource({{"KP", u0 * cox}}, in, 0.0,
                                          vdk, idk);
            checkTrue("U0-only and KP cards both converge", oku && okk);
            checkTrue("U0-only card conducts (U0 is no longer ignored)",
                      oku && std::abs(idu) > 1e-6);
            checkTrue("U0-only card matches the equivalent KP card",
                      oku && okk &&
                      std::abs(idu - idk) < 1e-6 * std::abs(idk));
        }

        // (c) Body effect: with K1 given the DC threshold must move with Vbs.
        {
            MosfetInstanceParams in;
            in.W = 2e-6; in.L = 0.5e-6;
            double v0 = 0.0, i0 = 0.0, vb = 0.0, ib = 0.0;
            const bool ok0 = commonSource({{"K1", 0.5}}, in, 0.0, v0, i0);
            const bool okb = commonSource({{"K1", 0.5}}, in, -1.0, vb, ib);
            checkTrue("body-bias sweep converges", ok0 && okb);
            (void)0;
            checkTrue("K1 produces a real body effect (Id falls at Vbs=-1V)",
                      ok0 && okb && ib < 0.9 * i0 && std::abs(ib) > 1e-9);
        }

        // (d) Junction depletion charge appears once CJS is non-zero.
        {
            const std::vector<KV> base = {{"VTH0", 1.5}, {"CJD", 0.0}, {"CJSWD", 0.0}, {"CJSWGD", 0.0}};
            MosfetInstanceParams in;
            in.W = 100e-6; in.L = 1e-6;
            in.AD = 1e-6; in.AS = 1e-6;
            double i0 = 0.0, q0 = 0.0, i1 = 0.0, q1 = 0.0;
            const bool ok0 = junctionOp(true, in, base, 1.0, i0, q0);
            std::vector<KV> withC = base;
            withC.push_back({"CJD", 5e-4});
            const bool ok1 = junctionOp(true, in, withC, 1.0, i1, q1);
            {
                ModelRegistry reg;
                RCP<const MosfetBSIM4Model> m0 = cardX(reg, "DA", true, base);
                ModelRegistry reg2;
                RCP<const MosfetBSIM4Model> m1 = cardX(reg2, "DB", true, withC);
                (void)m0; (void)m1;
            }
            checkTrue("CJS on/off both converge", ok0 && ok1);
            (void)0;
            checkTrue("junction depletion charge is zero with CJD=0",
                      ok0 && std::abs(q0) < 1e-13);
            checkTrue("junction depletion charge appears with CJD>0",
                      ok0 && ok1 && std::abs(q1 - q0) > 1e-15);
        }

        // (e) RDSMOD/RGATEMOD internal nodes are stamped, not floating.
        {
            MosfetInstanceParams in;
            in.W = 2e-6; in.L = 0.5e-6; in.NRD = 1.0; in.NRS = 1.0;
            double v0 = 0.0, i0 = 0.0, v1 = 0.0, i1 = 0.0;
            const bool ok0 = commonSource({{"RDSMOD", 1.0}}, in, 0.0, v0, i0);
            std::vector<KV> withR = {{"RDSMOD", 1.0}, {"RSH", 100.0}};
            const bool ok1 = commonSource(withR, in, 0.0, v1, i1);
            checkTrue("RDSMOD=1 converges with and without RSH",
                      ok0 && ok1);
            if (ok0 && ok1) {
                const double want = (3.0 - v1) / 1e3 - v1 / 1e6;
                checkTrue("KCL closes at the external drain node (RDSMOD=1)",
                          std::abs(want - i1) < 1e-9);
                checkTrue("series resistance reduces the drain current",
                          i1 <= i0 + 1e-12 && std::abs(i1) > 1e-10);
            }
            double vg0 = 0.0, ig0 = 0.0;
            const bool okg0 = commonSource({}, in, 0.0, vg0, ig0);
            double vg = 0.0, ig = 0.0;
            const bool okg = commonSource({{"RGATEMOD", 1.0}}, in, 0.0,
                                          vg, ig);
            checkTrue("RGATEMOD=1 converges (gate node not floating)", okg);
            if (okg && okg0) {
                checkTrue("RGATEMOD=1 leaves the drain current intact",
                          std::abs(ig - ig0) < 1e-3 * std::abs(ig0));
            }
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
