/**
 * @file test_mosfet.cpp
 * @brief Verification of the MOSFET level-1 (Shichman-Hodges) devices,
 *        NMOS + PMOS.
 *
 * Every check is analytic or a physics invariant:
 *  1. Registry mechanics (NMOS/PMOS, wrong-type + LEVEL=2 rejection).
 *  2. Threshold: off below VTO, on above.
 *  3. Saturation Id = (BETA/2)·(Vgs−Vth)²·(1+LAM·Vds).
 *  4. Linear Id = BETA·((Vgs−Vth)·Vds − Vds²/2)·(1+LAM·Vds).
 *  5. Transconductance gm vs finite difference.
 *  6. Body effect (GAMMA raises Vth, Id follows).
 *  7. Channel-length modulation (go = LAM·Id).
 *  8. PMOS mirror of NMOS.
 *  9. Reverse Vds (source/drain swap) matches the mirrored circuit.
 * 10. Bulk diode forward current.
 * 11. Cutoff leakage ~ 0.
 * 12. CMOS inverter transient (rail-to-rail).
 * 13. W scaling + sharing + temperature + parser errors.
 */

#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/core/NetlistParser.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/Capacitor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Mosfet.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/models/core_models/MosfetModel.h"
#include "../src/models/core/ModelRegistry.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"

#include <cmath>
#include <cstdio>
#include <string>

using namespace e_XSpice;

static int g_failures = 0;

static void check(const char* what, double got, double expect, double tol) {
    const double err = std::abs(got - expect);
    const bool ok = err <= tol;
    if (!ok) ++g_failures;
    std::printf("  [%s] %-40s got=% .8g expect=% .8g\n",
                ok ? "OK" : "FAIL", what, got, expect);
}

static void checkTrue(const char* what, bool cond) {
    if (!cond) ++g_failures;
    std::printf("  [%s] %s\n", cond ? "OK" : "FAIL", what);
}

static RCP<const MosfetModel> cardWith(ModelRegistry& reg,
                                       const std::string& name,
                                       const std::string& body, bool nmos) {
    ModelCard card;
    card.name = name;
    card.type = nmos ? "NMOS" : "PMOS";
    std::string token;
    for (char c : body + " ") {
        if (c == ' ' || c == '\t' || c == '\n') {
            if (!token.empty()) {
                const size_t eq = token.find('=');
                card.params[token.substr(0, eq)] =
                    std::stod(token.substr(eq + 1));
                token.clear();
            }
        } else {
            token.push_back(c);
        }
    }
    reg.createModel(card);
    return reg.findAs<MosfetModel>(name);
}

// BETA = KP·(W/L) with the test geometry below.
static constexpr double kW = 100e-6;
static constexpr double kL = 10e-6;

struct MosResult {
    bool ok = false;
    double id = 0, is = 0, ib = 0;
    double vth = 0, gm = 0;
    bool sat = false;
};

/**
 * @brief Common-source with forced terminal voltages: drain at vd,
 *        gate at vg, source at vs, bulk at vb (all absolute).
 */
static MosResult runForced(RCP<const MosfetModel> model, bool nmos,
                           double vd, double vg, double vs, double vb) {
    MosResult r;
    Circuit c("mos");
    Node* gnd = c.getGroundNode();
    const double p = nmos ? 1.0 : -1.0;

    // The PMOS mirror negates every supply (same N-oriented bias).
    Node* nd = c.addNode(1, "d");
    Node* ng = c.addNode(2, "g");
    Node* ns = c.addNode(3, "s");
    Node* nb = c.addNode(4, "b");
    c.addDevice(make_rcp<VoltageSource>("VD", nd, gnd,
                                        SourceWaveform::dc(p * vd)));
    c.addDevice(make_rcp<VoltageSource>("VG", ng, gnd,
                                        SourceWaveform::dc(p * vg)));
    c.addDevice(make_rcp<VoltageSource>("VS", ns, gnd,
                                        SourceWaveform::dc(p * vs)));
    c.addDevice(make_rcp<VoltageSource>("VB", nb, gnd,
                                        SourceWaveform::dc(p * vb)));
    MosfetInstanceParams inst;
    inst.W = kW;
    inst.L = kL;
    auto dev = make_rcp<Mosfet>("M1", nd, ng, ns, nb, model, inst);
    Mosfet* m = dev.get();
    c.addDevice(std::move(dev));

    DCAnalysis dc(c);
    DCResult res = dc.run(false, false);
    r.ok = res.success;
    if (!res.success) return r;
    r.id = m->getDrainCurrent();
    r.is = m->getSourceCurrent();
    r.ib = m->getBulkCurrent();
    r.vth = m->getVth();
    r.gm = m->getTransconductance();
    r.sat = m->inSaturation();
    return r;
}

// ===============================================================
// Tests
// ===============================================================

static void test_registry() {
    std::printf("\n1) Registry mechanics (NMOS/PMOS cards)\n");
    ModelRegistry reg;
    checkTrue("NMOS factory registered", reg.hasModelType("NMOS"));
    checkTrue("PMOS factory registered", reg.hasModelType("PMOS"));

    RCP<const MosfetModel> nm = cardWith(reg, "N1", "VTO=1 KP=2e-5", true);
    checkTrue("NMOS card created", (bool)nm);
    checkTrue("typeName NMOS", nm->typeName() == "NMOS");
    RCP<const MosfetModel> pm = cardWith(reg, "P1", "VTO=-1 KP=2e-5", false);
    checkTrue("PMOS card created", (bool)pm);
    checkTrue("typeName PMOS", pm->typeName() == "PMOS");

    bool badLevel = false;
    try {
        ModelCard c; c.name = "BL"; c.type = "NMOS"; c.params["LEVEL"] = 2.0;
        reg.createModel(c);
    } catch (const std::exception&) { badLevel = true; }
    checkTrue("LEVEL=2 rejected", badLevel);

    bool badKp = false;
    try {
        ModelCard c; c.name = "BK"; c.type = "NMOS"; c.params["KP"] = -1.0;
        reg.createModel(c);
    } catch (const std::exception&) { badKp = true; }
    checkTrue("negative KP rejected", badKp);

    {
        ModelCard c; c.name = "R1"; c.type = "R";
        reg.createModel(c);
        checkTrue("R card is not a MOSFET", !reg.findAs<MosfetModel>("R1"));
    }
}

static void test_threshold() {
    std::printf("\n2) Threshold behavior\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "TH", "VTO=1 KP=2e-5", true);

    const MosResult off = runForced(m, true, 5.0, 0.5, 0.0, 0.0);
    checkTrue("below Vth converges", off.ok);
    if (off.ok) {
        checkTrue("channel off (|Id| < 1nA)", std::abs(off.id) < 1e-9);
        checkTrue("not saturated", !off.sat);
    }
    const MosResult on = runForced(m, true, 5.0, 2.0, 0.0, 0.0);
    checkTrue("above Vth converges", on.ok);
    if (on.ok) {
        checkTrue("channel on (Id > 10uA)", on.id > 10e-6);
        checkTrue("saturated", on.sat);
    }
}

static void test_saturation() {
    std::printf("\n3) Saturation current\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "SAT", "VTO=1 KP=2e-5", true);

    // BETA = 20u·(100u/10u) = 200u; ov = 1 → Id = 100u·1 = 100uA.
    const MosResult r = runForced(m, true, 5.0, 2.0, 0.0, 0.0);
    checkTrue("DC converges", r.ok);
    if (!r.ok) return;
    check("Id = BETA/2·ov²", r.id, 0.5 * 200e-6 * 1.0, 1e-9);
    check("KCL: id+is+ib = 0", r.id + r.is + r.ib, 0.0, 1e-15);
    check("Vth = VTO (no body)", r.vth, 1.0, 1e-12);
}

static void test_linear() {
    std::printf("\n4) Linear-region current\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "LIN", "VTO=1 KP=2e-5", true);

    // ov = 1, vds = 0.2 → Id = 200u·(0.2 − 0.02) = 36uA.
    const MosResult r = runForced(m, true, 0.2, 2.0, 0.0, 0.0);
    checkTrue("DC converges", r.ok);
    if (!r.ok) return;
    check("Id linear formula", r.id, 200e-6 * (0.2 - 0.02), 1e-9);
    checkTrue("not saturated", !r.sat);
}

static void test_gm() {
    std::printf("\n5) Transconductance vs finite difference\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "GM", "VTO=1 KP=2e-5", true);

    const MosResult r1 = runForced(m, true, 5.0, 2.0, 0.0, 0.0);
    const MosResult r2 = runForced(m, true, 5.0, 2.01, 0.0, 0.0);
    checkTrue("both converge", r1.ok && r2.ok);
    if (!r1.ok || !r2.ok) return;
    const double gmFd = (r2.id - r1.id) / 0.01;
    check("gm = BETA·ov", r1.gm, 200e-6 * 1.0, 1e-9);
    check("gm matches Id slope", gmFd, r1.gm, 0.01 * r1.gm);
}

static void test_body_effect() {
    std::printf("\n6) Body effect\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(
        reg, "BODY", "VTO=1 KP=2e-5 GAMMA=0.5 PHI=0.6", true);

    // Bulk at −2 V (reverse): Vth = 1 + 0.5·(√2.6 − √0.6).
    const MosResult r = runForced(m, true, 5.0, 2.0, 0.0, -2.0);
    checkTrue("DC converges", r.ok);
    if (!r.ok) return;
    const double vthExpect =
        1.0 + 0.5 * (std::sqrt(2.6) - std::sqrt(0.6));
    check("Vth with body effect", r.vth, vthExpect, 1e-9);
    const double ov = 2.0 - vthExpect;
    check("Id follows Vth", r.id, 0.5 * 200e-6 * ov * ov, 2e-9);
}

static void test_lambda() {
    std::printf("\n7) Channel-length modulation\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(
        reg, "LAM", "VTO=1 KP=2e-5 LAMBDA=0.02", true);

    const MosResult r3 = runForced(m, true, 3.0, 2.0, 0.0, 0.0);
    const MosResult r5 = runForced(m, true, 5.0, 2.0, 0.0, 0.0);
    checkTrue("both converge", r3.ok && r5.ok);
    if (!r3.ok || !r5.ok) return;
    const double go = (r5.id - r3.id) / 2.0;
    checkTrue("go > 0", go > 0.0);
    // Id = Id0·(1+LAM·Vds) ⇒ go = LAM·Id0 is Vds-independent.
    check("go = LAM·Id0", go, 0.02 * r3.id / (1.0 + 0.02 * 3.0), 1e-12);
}

static void test_pmos_mirror() {
    std::printf("\n8) PMOS mirror of NMOS\n");
    ModelRegistry reg;
    RCP<const MosfetModel> nm = cardWith(reg, "NM", "VTO=1 KP=2e-5 LAMBDA=0.02",
                                         true);
    RCP<const MosfetModel> pm = cardWith(reg, "PM", "VTO=-1 KP=2e-5 LAMBDA=0.02",
                                         false);

    const MosResult n = runForced(nm, true, 5.0, 2.0, 0.0, 0.0);
    const MosResult p = runForced(pm, false, 5.0, 2.0, 0.0, 0.0);
    checkTrue("both converge", n.ok && p.ok);
    if (!n.ok || !p.ok) return;
    check("PMOS id = −NMOS id", p.id, -n.id, 1e-12);
    check("PMOS is = −NMOS is", p.is, -n.is, 1e-12);
    check("PMOS ib = −NMOS ib", p.ib, -n.ib, 1e-15);
    check("PMOS Vth(n) = NMOS Vth", p.vth, n.vth, 1e-12);
}

static void test_reverse_vds() {
    std::printf("\n9) Reverse Vds (source/drain swap)\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "REV", "VTO=1 KP=2e-5", true);

    // Physical: drain 4.8 V, source 5 V (vds = −0.2), gate 7 V, bulk 5 V.
    // Effective: vgs_e = 2.2 − ... vgd = 2.2? vgs=2, vds=−0.2 → vgs_e = 2.2,
    // vds_e = 0.2, ov = 1.2 → Id_eff = 200u·(1.2·0.2 − 0.02) = 44uA into
    // the physical SOURCE, i.e. physical drain current −44uA.
    const MosResult r = runForced(m, true, 4.8, 7.0, 5.0, 5.0);
    checkTrue("DC converges", r.ok);
    if (!r.ok) return;
    check("swapped Id", r.id, -200e-6 * (1.2 * 0.2 - 0.02), 2e-9);
    check("KCL", r.id + r.is + r.ib, 0.0, 1e-15);

    // The mirrored forward circuit (supplies exchanged) must agree.
    const MosResult f = runForced(m, true, 5.0, 7.0, 4.8, 4.8);
    checkTrue("forward twin converges", f.ok);
    if (f.ok) {
        check("swap symmetry: id_rev = −id_fwd", r.id, -f.id, 2e-9);
    }
}

static void test_bulk_diode() {
    std::printf("\n10) Bulk diode forward + 11) cutoff leakage\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "BD", "VTO=1 KP=2e-5 IS=1e-14",
                                        true);

    // Bulk 0.7 V above source (forward) with the drain high so only the
    // B–S diode is forward (B–D stays reverse). Tolerance covers the
    // Xyce-vs-CODATA physical-constant difference in Vt (2e-5 relative
    // in k/q amplifies through exp(27)).
    const MosResult f = runForced(m, true, 5.0, 0.0, 0.0, 0.7);
    checkTrue("forward-bulk converges", f.ok);
    if (f.ok) {
        const double vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
        const double expect = 1e-14 * (std::exp(0.7 / vt) - 1.0);
        check("Ibs forward", -f.ib, expect, 1e-2 * expect);
    }

    const MosResult c = runForced(m, true, 5.0, 0.0, 0.0, 0.0);
    checkTrue("cutoff converges", c.ok);
    if (c.ok) {
        checkTrue("cutoff leakage < 1nA", std::abs(c.id) < 1e-9);
    }
}

static void test_cmos_inverter() {
    std::printf("\n12) CMOS inverter: DC transfer + transient (rail-to-rail)\n");
    ModelRegistry reg;
    // Realistic overlap + bulk caps: they load the output (companion anchor
    // for Newton) and exercise the Meyer/overlap/bulk charge code. A
    // capacitance-free output node admits a second, spurious DC solution
    // (PMOS bulk diode balancing the NMOS drain current at VDD+0.7); the
    // DC transfer below pins the physical branch with C=0 cards instead.
    RCP<const MosfetModel> nm = cardWith(
        reg, "NI", "VTO=1 KP=2e-5 CGSO=4e-10 CGDO=4e-10 CBS=20e-15 CBD=20e-15",
        true);
    RCP<const MosfetModel> pm = cardWith(
        reg, "PI", "VTO=-1 KP=2e-5 CGSO=4e-10 CGDO=4e-10 CBS=20e-15 CBD=20e-15",
        false);

    // DC transfer with capacitance-free cards, swept upward WITH warm
    // start: each point starts from its converged neighbor, so Newton
    // tracks the physical branch (cold starts may land on the parasitic
    // bulk-diode branch — a genuine second solution, not a solver bug).
    {
        Circuit c("invdc");
        Node* gnd = c.getGroundNode();
        Node* nvdd = c.addNode(1, "vdd");
        Node* nin = c.addNode(2, "in");
        Node* nout = c.addNode(3, "out");
        c.addDevice(make_rcp<VoltageSource>("VDD", nvdd, gnd,
                                            SourceWaveform::dc(5.0)));
        auto vinSrc = make_rcp<VoltageSource>("VIN", nin, gnd,
                                              SourceWaveform::dc(0.0));
        VoltageSource* vinp = vinSrc.get();
        c.addDevice(std::move(vinSrc));
        ModelRegistry reg2;
        RCP<const MosfetModel> nm0 = cardWith(reg2, "N0", "VTO=1 KP=2e-5", true);
        RCP<const MosfetModel> pm0 = cardWith(reg2, "P0", "VTO=-1 KP=2e-5",
                                              false);
        MosfetInstanceParams inst;
        inst.W = kW;
        inst.L = kL;
        c.addDevice(make_rcp<Mosfet>("MN", nout, nin, gnd, gnd, nm0, inst));
        c.addDevice(make_rcp<Mosfet>("MP", nout, nin, nvdd, nvdd, pm0, inst));

        DCAnalysis dc(c);
        auto dcAt = [&](double vin) {
            vinp->setWaveform(SourceWaveform::dc(vin));
            DCResult r = dc.run(false, true);   // warm-started continuation
            return r.success ? nout->getVoltage() : 1e99;
        };
        check("DC transfer: in=0 → out=VDD", dcAt(0.0), 5.0, 1e-4);
        // High-gain region (LAM=0 ⇒ vertical transfer): both devices
        // saturated over out ∈ [vin−1, vin+1]; any point there is valid.
        double vmid = 1e99;
        for (double vin = 0.25; vin <= 2.51; vin += 0.25) vmid = dcAt(vin);
        checkTrue("DC transfer: mid-point in saturated window",
                  vmid > 1.0 && vmid < 4.0);
        double vout = 1e99;
        for (double vin = 2.75; vin <= 5.01; vin += 0.25) vout = dcAt(vin);
        check("DC transfer: in=5 → out=0", vout, 0.0, 0.15);
    }

    Circuit c("inv");
    Node* gnd = c.getGroundNode();
    Node* nvdd = c.addNode(1, "vdd");
    Node* nin = c.addNode(2, "in");
    Node* nout = c.addNode(3, "out");
    c.addDevice(make_rcp<VoltageSource>("VDD", nvdd, gnd,
                                        SourceWaveform::dc(5.0)));
    c.addDevice(make_rcp<VoltageSource>(
        "VIN", nin, gnd,
        SourceWaveform::pulse(0.0, 5.0, 1e-6, 1e-8, 1e-8, 2e-6, 4e-6)));
    MosfetInstanceParams inst;
    inst.W = kW;
    inst.L = kL;
    c.addDevice(make_rcp<Mosfet>("MN", nout, nin, gnd, gnd, nm, inst));
    c.addDevice(make_rcp<Mosfet>("MP", nout, nin, nvdd, nvdd, pm, inst));
    // Realistic load: slows the output (C·V/I ≈ 3 ns) so the 10 ns input
    // edge is resolved instead of jumped over in a single step.
    c.addDevice(make_rcp<Capacitor>("CL", nout, gnd, 1e-12));

    const auto& nodes = c.getNodes();
    auto idx = [&](Node* n) {
        for (size_t i = 0; i < nodes.size(); ++i)
            if (nodes[i].get() == n) return static_cast<int>(i);
        return -1;
    };
    const int iout = idx(nout);

    TransientAnalysis tran(c);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    TransientResult r = tran.run(0.0, 4e-6, 2e-9, false);
    checkTrue("inverter transient converges", r.success);
    if (!r.success) {
        std::printf("  [diag] %s\n", r.message.c_str());
        return;
    }
    if (r.data.empty()) return;

    // t = 0: in = 0 → out = VDD. End (t = 4us): back to 0 → out = VDD.
    check("out high at t=0", r.data.front().node_voltages[iout], 5.0, 0.05);
    // Mid-high input (t ≈ 2us): out = 0.
    double vmid = 5.0;
    double tmid = 0.0;
    int nwin = 0;
    for (const auto& p : r.data) {
        if (p.time >= 1.5e-6 && p.time <= 2.5e-6) {
            ++nwin;
            if (p.node_voltages[iout] < vmid) {
                vmid = p.node_voltages[iout];
                tmid = p.time;
            }
        }
    }
    std::printf("  (points=%zu in-window=%d lowest out: %.6g V at t=%.6g)\n",
                r.data.size(), nwin, vmid, tmid);
    checkTrue("out pulled low during input high", vmid < 0.1);
}

static void test_scaling_sharing_temp_parser() {
    std::printf("\n13) W scaling, sharing, temperature, parser\n");
    ModelRegistry reg;
    RCP<const MosfetModel> m = cardWith(reg, "SC", "VTO=1 KP=2e-5", true);

    // Double width → double current.
    auto runW = [&](double w) {
        Circuit c("w");
        Node* gnd = c.getGroundNode();
        Node* nd = c.addNode(1, "d");
        Node* ng = c.addNode(2, "g");
        c.addDevice(make_rcp<VoltageSource>("VD", nd, gnd,
                                            SourceWaveform::dc(5.0)));
        c.addDevice(make_rcp<VoltageSource>("VG", ng, gnd,
                                            SourceWaveform::dc(2.0)));
        MosfetInstanceParams inst;
        inst.W = w;
        inst.L = kL;
        auto dev = make_rcp<Mosfet>("M1", nd, ng, gnd, gnd, m, inst);
        Mosfet* mp = dev.get();
        c.addDevice(std::move(dev));
        DCAnalysis dc(c);
        DCResult r = dc.run(false, false);
        return r.success ? mp->getDrainCurrent() : -1.0;
    };
    const double i1 = runW(kW);
    const double i2 = runW(2.0 * kW);
    checkTrue("both converge", i1 > 0.0 && i2 > 0.0);
    if (i1 > 0.0 && i2 > 0.0) {
        check("Id(2W)/Id(W) = 2", i2 / i1, 2.0, 1e-6);
    }

    // Temperature: hotter card ⇒ larger bulk IS(T).
    RCP<const MosfetModel> hot = cardWith(reg, "HT", "VTO=1 KP=2e-5 TEMP=125",
                                         true);
    checkTrue("hot IS(T) > cold IS(T)",
              hot->thermal().isat > m->thermal().isat);

    // Netlist end-to-end.
    {
        Circuit c("net");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL NM NMOS (VTO=1 KP=2e-5)\n"
            ".MODEL PM PMOS (VTO=-1 KP=2e-5)\n"
            "VDD 1 0 DC 5\n"
            "VIN 2 0 DC 2\n"
            "M1 3 2 0 0 NM W=100u L=10u\n"
            "M2 3 2 1 1 PM W=100u L=10u\n"
            "RD 1 3 1k\n"
            ".OP\n.END\n",
            c);
        checkTrue("MOS netlist parses", ok);
        if (ok) {
            Mosfet* m1 = static_cast<Mosfet*>(c.getDevice("M1"));
            check("W parsed", m1->instanceParams().W, 100e-6, 1e-12);
            check("L parsed", m1->instanceParams().L, 10e-6, 1e-16);
            DCAnalysis dc(c);
            DCResult r = dc.run(false, false);
            checkTrue("MOS netlist DC converges", r.success);
        }
    }
    {
        Circuit c("e1");
        NetlistParser p;
        const bool ok = p.parseString("M1 1 2 0 MM\n.END\n", c);
        checkTrue("3-terminal M rejected", !ok);
    }
    {
        Circuit c("e2");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-14)\nM1 1 2 0 0 DM\n.END\n", c);
        checkTrue("diode card rejected for M1", !ok);
    }
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("=== MOSFET level-1 (Shichman-Hodges) verification ===\n");
    test_registry();
    test_threshold();
    test_saturation();
    test_linear();
    test_gm();
    test_body_effect();
    test_lambda();
    test_pmos_mirror();
    test_reverse_vds();
    test_bulk_diode();
    test_cmos_inverter();
    test_scaling_sharing_temp_parser();

    std::printf("\n=== %s (%d failure%s) ===\n",
                g_failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    mpi.finalize();
    return g_failures == 0 ? 0 : 1;
}
