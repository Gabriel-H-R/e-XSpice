/**
 * @file test_bjt.cpp
 * @brief Verification of the BJT devices (NPN/PNP) and their model card.
 *
 * Every check is either analytic (closed form / reference bisection) or a
 * physics invariant, so a regression in the Ebers-Moll laws, the charge
 * stamping or the polarity handling is caught here:
 *
 *  1. Registry mechanics for NPN/PNP cards (incl. wrong-type rejection).
 *  2. Card aliases (VA, VB, CJ, PE, ME, PC, MC; CJS is substrate-only).
 *  3. Active region:  i_c/i_b = BF exactly, Vce, Vbe vs bisection.
 *  4. Early effect:  ic(Vce) ratio + output conductance sign.
 *  5. Saturation:    Vce_sat small, ic < BF·ib.
 *  6. Cutoff:        collector pulled to the rail (gmin-limited).
 *  7. PNP:           exact mirror of the NPN solution (signs flip).
 *  8. AREA:          currents scale, betas do not change.
 *  9. Temperature:   IS(T) up => Vbe down (fixed base current).
 * 10. Charge model:  Q_BE/Q_BC equal the closed-form card value.
 * 11. Transient:     cutoff -> active -> saturation endpoints.
 * 12. Sharing:       two transistors on one card.
 * 13. Parser:        netlist end-to-end + error paths.
 */

#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/core/NetlistParser.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/CurrentSource.h"
#include "../src/devices/Bjt.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/models/core_models/BjtModel.h"
#include "../src/models/core/ModelRegistry.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"

#include <cmath>
#include <cstdio>
#include <string>

using namespace e_XSpice;

static int g_failures = 0;
static const double kQ = 1.602176634e-19;
static const double kK = 1.380649e-23;
static const double kVt = kK * 300.15 / kQ;   // thermal voltage at TNOM

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

/**
 * @brief Builds a card through the registry from a "KEY=value ..." body.
 */
static RCP<const BjtModel> cardWith(ModelRegistry& reg, const std::string& name,
                                   const std::string& body, bool pnp) {
    ModelCard card;
    card.name = name;
    card.type = pnp ? "PNP" : "NPN";
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
    return reg.findAs<BjtModel>(name);
}

// ---------------------------------------------------------------
// Reference circuits
// ---------------------------------------------------------------

/**
 * @brief Common-emitter bias: VCC - RC - collector, forced base current.
 *        Returns false when the analysis fails.
 */
struct CeResult {
    bool ok = false;
    double vbe = 0, vbc = 0, vce = 0, vc = 0;
    double ib = 0, ic = 0, ie = 0;
    double qbe = 0, qbc = 0;   // charges (copied out while the circuit lives)
};

static CeResult runCe(RCP<const BjtModel> model, bool pnp, double vcc,
                      double rc, double ibias, double area = 1.0) {
    CeResult r;
    r.ok = false;
    Circuit c("ce");
    Node* gnd = c.getGroundNode();
    Node* nrail = c.addNode(1, "1");
    Node* ncol = c.addNode(2, "2");
    Node* nbase = c.addNode(3, "3");

    // PNP mirror: everything negated (rail at -vcc, base current drawn).
    const double p = pnp ? -1.0 : 1.0;
    c.addDevice(make_rcp<VoltageSource>("V1", nrail, gnd,
                                        SourceWaveform::dc(p * vcc)));
    c.addDevice(make_rcp<Resistor>("RC", nrail, ncol, rc));
    c.addDevice(make_rcp<CurrentSource>("IB", p > 0 ? gnd : nbase,
                                        p > 0 ? nbase : gnd,
                                        SourceWaveform::dc(ibias)));
    BjtInstanceParams inst;
    inst.AREA = area;
    c.addDevice(make_rcp<Bjt>("Q1", ncol, nbase, gnd, model, inst));

    DCAnalysis dc(c);
    DCResult res = dc.run(false, false);
    r.ok = res.success;
    if (!res.success) {
        std::printf("  [diag] DC failed: %s\n", res.message.c_str());
        return r;
    }
    Bjt* dev = static_cast<Bjt*>(c.getDevice("Q1"));
    r.vc = ncol->getVoltage();
    r.vbe = nbase->getVoltage();
    r.vbc = nbase->getVoltage() - r.vc;
    r.vce = r.vc;
    r.ib = dev->getBaseCurrent();
    r.ic = dev->getCollectorCurrent();
    r.ie = dev->getEmitterCurrent();
    r.qbe = dev->getBaseCharge();
    r.qbc = dev->getCollectorCharge();
    return r;
}

/**
 * @brief Reference base current of the card law at vbe:
 *        i_b = beta_F·(exp(vbe/vte) - 1) + gmin·vbe
 */
static double baseCurrentAt(const BjtModel& m, double vbe) {
    const double beIsat = m.thermal().beIsat;
    const double betaF = beIsat / m.BF;
    return betaF * (std::exp(vbe / m.thermal().vte) - 1.0) + m.GMIN * vbe;
}

// ===============================================================
// Tests
// ===============================================================

static void test_registry() {
    std::printf("\n1) Registry mechanics (NPN/PNP cards)\n");
    ModelRegistry reg;
    checkTrue("NPN factory registered", reg.hasModelType("NPN"));
    checkTrue("PNP factory registered", reg.hasModelType("PNP"));
    checkTrue("D factory still registered", reg.hasModelType("D"));

    RCP<const BjtModel> npn = cardWith(reg, "QN", "IS=1e-16 BF=100", false);
    checkTrue("NPN card created", (bool)npn);
    checkTrue("typeName NPN", npn->typeName() == "NPN");
    checkTrue("not pnp", !npn->isPnp());

    RCP<const BjtModel> pnp = cardWith(reg, "QP", "IS=1e-16 BF=100", true);
    checkTrue("PNP card created", (bool)pnp);
    checkTrue("typeName PNP", pnp->typeName() == "PNP");
    checkTrue("is pnp", pnp->isPnp());

    bool dup = false;
    try { cardWith(reg, "QN", "IS=1e-16", false); }
    catch (const std::exception&) { dup = true; }
    checkTrue("duplicate name rejected", dup);

    bool badBeta = false;
    try {
        ModelCard card; card.name = "BAD"; card.type = "NPN";
        card.params["BF"] = -1.0;
        reg.createModel(card);
    } catch (const std::exception&) { badBeta = true; }
    checkTrue("invalid BF rejected", badBeta);

    bool badType = false;
    try {
        ModelCard card; card.name = "X"; card.type = "Q";
        reg.createModel(card);
    } catch (const std::exception&) { badType = true; }
    checkTrue("unknown type Q rejected", badType);

    // A resistor card must NOT be retrievable as a BJT.
    {
        ModelCard card; card.name = "R1"; card.type = "R";
        reg.createModel(card);
        checkTrue("R card is not a BJT", !reg.findAs<BjtModel>("R1"));
    }
}

static void test_aliases() {
    std::printf("\n2) Card aliases\n");
    ModelRegistry reg;
    RCP<const BjtModel> m = cardWith(
        reg, "AL",
        "VA=50 VB=200 CJ=5e-12 PE=0.7 ME=0.4 CJC=3e-12 PC=0.6 MC=0.35", false);
    check("VAF via VA", m->VAF, 50.0, 0.0);
    check("VAR via VB", m->VAR, 200.0, 0.0);
    check("CJE via CJ", m->CJE, 5e-12, 0.0);
    check("VJE via PE", m->VJE, 0.7, 0.0);
    check("MJE via ME", m->MJE, 0.4, 0.0);
    check("CJC direct", m->CJC, 3e-12, 0.0);
    check("VJC via PC", m->VJC, 0.6, 0.0);
    check("MJC via MC", m->MJC, 0.35, 0.0);
    checkTrue("thermal scaled", m->thermal().vt > 0.02);
    checkTrue("beIsat > 0", m->thermal().beIsat > 0.0);
}

static void test_active_region() {
    std::printf("\n3) Active region: beta, Vce, Vbe\n");
    ModelRegistry reg;
    RCP<const BjtModel> m =
        cardWith(reg, "ACT", "IS=1e-16 BF=100 NF=1", false);

    const CeResult r = runCe(m, false, 10.0, 1000.0, 10e-6);
    checkTrue("DC converges", r.ok);
    if (!r.ok) return;

    check("beta = ic/ib", r.ic / r.ib, 100.0, 1e-3);
    check("KCL: ic = -(ib+ie)", r.ib + r.ie + r.ic, 0.0, 1e-15);

    // Collector voltage: v_c = VCC - ic·RC (i_c>0 flows into the rail side).
    check("Vce = VCC - ic*RC", r.vc, 10.0 - r.ic * 1000.0, 1e-5);

    // Vbe from the reference bisection of the base-current law.
    double lo = 0.0, hi = 1.5;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (baseCurrentAt(*m, mid) < 10e-6) lo = mid; else hi = mid;
    }
    check("Vbe vs reference", r.vbe, 0.5 * (lo + hi), 1e-7);
}

static void test_early() {
    std::printf("\n4) Early effect (base-width modulation)\n");
    ModelRegistry reg;
    RCP<const BjtModel> m =
        cardWith(reg, "EAR", "IS=1e-16 BF=100 VAF=100", false);

    const CeResult r5 = runCe(m, false, 5.0, 1000.0, 10e-6);
    const CeResult r10 = runCe(m, false, 10.0, 1000.0, 10e-6);
    checkTrue("both converges", r5.ok && r10.ok);
    if (!r5.ok || !r10.ok) return;

    // Fixed base current => Vbe identical => ic ratio is purely Early:
    //   ic(Vce) ∝ 1 + (Vce - Vbe)/VAF
    check("Vbe identical at both rails", r10.vbe, r5.vbe, 1e-5);
    const double k5 = 1.0 + (r5.vce - r5.vbe) / 100.0;
    const double k10 = 1.0 + (r10.vce - r10.vbe) / 100.0;
    check("ic ratio follows Early model", r10.ic / r5.ic, k10 / k5, 1e-3);

    // Output conductance: the model is linear in Vce, so the finite
    // difference over the interval is exactly IS·exp(vbe/vte)/VAF = ic/k.
    const double go = (r10.ic - r5.ic) / (r10.vce - r5.vce);
    const double goExpect = r5.ic / (100.0 * k5);
    check("output conductance ~ ic/(VAF*k)", go, goExpect, 0.01 * goExpect);
    checkTrue("g_o > 0", go > 0.0);
}

static void test_saturation_and_cutoff() {
    std::printf("\n5) Saturation\n");
    ModelRegistry reg;
    RCP<const BjtModel> m =
        cardWith(reg, "SAT", "IS=1e-16 BF=100", false);

    // 1mA base current, 1k load, 5V rail: the load could take 5mA, the
    // transistor BF*ib = 100mA => deep saturation.
    const CeResult s = runCe(m, false, 5.0, 1000.0, 1e-3);
    checkTrue("saturation converges", s.ok);
    if (!s.ok) return;
    checkTrue("Vce_sat < 0.3V", s.vce < 0.3 && s.vce > 0.0);
    checkTrue("ic < BF*ib (not active)", s.ic < 100.0 * s.ib);
    checkTrue("both junctions forward", s.vbe > 0.0 && s.vbc > 0.0);

    std::printf("\n6) Cutoff\n");
    const CeResult c = runCe(m, false, 5.0, 1000.0, -1e-6);
    checkTrue("cutoff converges", c.ok);
    if (!c.ok) return;
    // Base pulled backwards: both junctions off, only gmin leaks. The rail
    // current is below 1uA, so v_c stays within millivolts of 5V.
    check("Vce ~ VCC in cutoff", c.vce, 5.0, 5e-3);
    checkTrue("collector current ~ 0", std::abs(c.ic) < 1e-6);
    checkTrue("emitter KCL", std::abs(c.ib + c.ie + c.ic) < 1e-15);
}

static void test_pnp_mirror() {
    std::printf("\n7) PNP mirror of the NPN solution\n");
    ModelRegistry reg;
    RCP<const BjtModel> npn = cardWith(reg, "QN", "IS=1e-16 BF=100", false);
    RCP<const BjtModel> pnp = cardWith(reg, "QP", "IS=1e-16 BF=100", true);

    const CeResult n = runCe(npn, false, 10.0, 1000.0, 10e-6);
    const CeResult p = runCe(pnp, true, 10.0, 1000.0, 10e-6);
    checkTrue("both converge", n.ok && p.ok);
    if (!n.ok || !p.ok) return;

    check("PNP Vbe = -NPN Vbe", p.vbe, -n.vbe, 1e-7);
    check("PNP Vbc = -NPN Vbc", p.vbc, -n.vbc, 1e-7);
    check("PNP Vc  = -NPN Vc", p.vc, -n.vc, 1e-7);
    check("PNP ic  = -NPN ic", p.ic, -n.ic, 1e-12);
    check("PNP ib  = -NPN ib", p.ib, -n.ib, 1e-15);
    check("PNP ie = -NPN ie", p.ie, -n.ie, 1e-12);

    // PNP saturation too (mirrored): deep saturation with a large base pull.
    const CeResult ps = runCe(pnp, true, 5.0, 1000.0, 1e-3);
    checkTrue("PNP saturation converges", ps.ok);
    if (ps.ok) {
        checkTrue("PNP Vce_sat small", ps.vce > -0.3 && ps.vce < 0.0);
        checkTrue("PNP Vbe < 0", ps.vbe < 0.0);
    }
}

static void test_area() {
    std::printf("\n8) AREA scaling\n");
    ModelRegistry reg;
    RCP<const BjtModel> m =
        cardWith(reg, "AR", "IS=1e-16 BF=100 CJC=1p", false);

    const CeResult a1 = runCe(m, false, 10.0, 1000.0, 10e-6, 1.0);
    const CeResult a2 = runCe(m, false, 10.0, 1000.0, 10e-6, 2.0);
    checkTrue("both converge", a1.ok && a2.ok);
    if (!a1.ok || !a2.ok) return;

    // The forced base current is fixed, so Vbe must SHIFT so that the
    // doubled saturation current still carries the same current:
    //   beta_F = IS(T)/BF  ->  vbe(2·area) = vbe(1) − vte·ln(2)
    check("Vbe shift = -vte*ln2", a2.vbe - a1.vbe,
          -m->thermal().vte * std::log(2.0), 1e-5);
    check("beta unchanged", a2.ic / a2.ib, a1.ic / a1.ib, 1e-3);
}

static void test_temperature() {
    std::printf("\n9) Temperature propagation\n");
    ModelRegistry reg;
    RCP<const BjtModel> cold = cardWith(reg, "TC", "IS=1e-16 BF=100", false);
    RCP<const BjtModel> hot = cardWith(reg, "TH", "IS=1e-16 BF=100 TEMP=125", false);

    checkTrue("hot IS(T) > cold IS(T)",
              hot->thermal().beIsat > cold->thermal().beIsat);

    const CeResult c = runCe(cold, false, 10.0, 1000.0, 10e-6);
    const CeResult h = runCe(hot, false, 10.0, 1000.0, 10e-6);
    checkTrue("both converge", c.ok && h.ok);
    if (!c.ok || !h.ok) return;
    // Hotter card => larger IS => smaller forward drop at the same current.
    checkTrue("Vbe decreases with T", h.vbe < c.vbe);
    std::printf("  (cold=%.6f V hot=%.6f V)\n", c.vbe, h.vbe);
}

static void test_charge_model() {
    std::printf("\n10) Charge model (Q_BE / Q_BC closed form)\n");
    ModelRegistry reg;
    RCP<const BjtModel> m = cardWith(
        reg, "CH", "IS=1e-16 BF=100 TF=1e-9 CJC=2e-12 VJC=0.75 MJC=0.33",
        false);

    const CeResult r = runCe(m, false, 10.0, 1000.0, 10e-6);
    checkTrue("DC converges", r.ok);
    if (!r.ok) return;

    // Independent closed-form evaluation of the card's charge partition.
    const double area = 1.0;
    const JunctionGeom geomBC{m->thermal().bcCap.cjo * area, 0.0, m->MJC};
    const double qdepBC = junctionCharge(r.vbc, 0.0, 0.0, geomBC,
                                         m->thermal().bcCap).q;
    check("Q_BC = dep(vbc) - TR*ic (TR=0)", r.qbc,
          qdepBC, 1e-7 * std::abs(qdepBC) + 1e-21);

    const JunctionGeom geomBE{m->thermal().beCap.cjo * area, 0.0, m->MJE};
    const double qdepBE = junctionCharge(r.vbe, 0.0, 0.0, geomBE,
                                         m->thermal().beCap).q;
    const double qbeExpect = qdepBE + m->TF * r.ic;
    check("Q_BE = dep(vbe) + TF*ic", r.qbe,
          qbeExpect, 1e-7 * std::abs(qbeExpect) + 1e-21);

    // Sanity: with TF=1ns and ic=1mA the transit charge dominates.
    checkTrue("transit charge present",
              std::abs(r.qbe) > 0.8e-12);
}

static void test_transient_switching() {
    std::printf("\n11) Transient switching (cutoff -> active)\n");
    ModelRegistry reg;
    RCP<const BjtModel> m =
        cardWith(reg, "SW", "IS=1e-16 BF=100 CJC=1e-12 CJE=2e-12", false);

    Circuit c("sw");
    Node* gnd = c.getGroundNode();
    Node* nrail = c.addNode(1, "1");      // 10 V rail
    Node* ndrv = c.addNode(3, "3");       // base drive source
    Node* nbase = c.addNode(4, "4");      // base terminal
    Node* ncol = c.addNode(2, "2");       // collector

    c.addDevice(make_rcp<VoltageSource>("VCC", nrail, gnd,
                                        SourceWaveform::dc(10.0)));
    c.addDevice(make_rcp<Resistor>("RC", nrail, ncol, 1000.0));
    c.addDevice(make_rcp<VoltageSource>(
        "VB", ndrv, gnd,
        SourceWaveform::pulse(0.0, 5.0, 1e-6, 1e-8, 1e-8, 1e-3, 1e-3)));
    c.addDevice(make_rcp<Resistor>("RB", ndrv, nbase, 100000.0));
    c.addDevice(make_rcp<Bjt>("Q1", ncol, nbase, gnd, m));

    // node_voltages are indexed by position in getNodes(), not by node id.
    const auto& nodes = c.getNodes();
    auto idx = [&](Node* n) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (nodes[i].get() == n) return static_cast<int>(i);
        }
        return -1;
    };
    const int icol = idx(ncol);
    const int ibase = idx(nbase);

    TransientAnalysis tran(c);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    const double dt = 2e-8;
    TransientResult r = tran.run(0.0, 2e-6, dt, false);
    checkTrue("transient converges", r.success);
    if (!r.success || r.data.empty()) return;

    // t = 0: base at 0 V => cutoff, the collector sits on the rail.
    check("cutoff: Vc = VCC", r.data.front().node_voltages[icol], 10.0, 5e-3);

    // End of the pulse: base at 5 V through 100k => active region with
    // ic = BF·ib = BF·(5 - Vbe)/RB.
    const auto& last = r.data.back();
    const double vbe = last.node_voltages[ibase];
    const double ib = (5.0 - vbe) / 100000.0;
    const double icExpect = 100.0 * ib;
    const double vc = last.node_voltages[icol];
    check("final ic = BF*ib", (10.0 - vc) / 1000.0, icExpect,
          0.02 * icExpect);
    checkTrue("switched to the active region (Vce < VCC - 1)", vc < 9.0);

    // Turn-on: the base charges through RB (100k) with tau ~ 200 ns, so the
    // collector is still on the rail just after the edge and is down by the
    // end of the run. Nothing may swing above its cutoff value.
    double vcEarly = r.data.front().node_voltages[icol];
    for (const auto& p : r.data) {
        if (p.time > 1.1e-6) break;
        vcEarly = p.node_voltages[icol];
    }
    checkTrue("turn-on started by 1.1us", vcEarly < 9.0);
    checkTrue("collector switched on by the end",
              last.node_voltages[icol] < vcEarly - 0.5);
    std::printf("  (final Vbe=%.6f V, Vce=%.6f V, ib=%.6g A)\n",
                vbe, vc, ib);
}

static void test_sharing_and_parser() {
    std::printf("\n12) Sharing + 13) parser\n");
    ModelRegistry reg;
    RCP<const BjtModel> m = cardWith(reg, "SH", "IS=1e-16 BF=100", false);

    {
        Circuit c("sh");
        Node* gnd = c.getGroundNode();
        Node* nr1 = c.addNode(1, "1");
        Node* nc1 = c.addNode(2, "2");
        Node* nb1 = c.addNode(3, "3");
        Node* nr2 = c.addNode(11, "11");
        Node* nc2 = c.addNode(12, "12");
        Node* nb2 = c.addNode(13, "13");
        c.addDevice(make_rcp<VoltageSource>("V1", nr1, gnd,
                                            SourceWaveform::dc(10.0)));
        c.addDevice(make_rcp<VoltageSource>("V2", nr2, gnd,
                                            SourceWaveform::dc(10.0)));
        c.addDevice(make_rcp<Resistor>("RC1", nr1, nc1, 1000.0));
        c.addDevice(make_rcp<Resistor>("RC2", nr2, nc2, 1000.0));
        c.addDevice(make_rcp<CurrentSource>("IB1", gnd, nb1,
                                            SourceWaveform::dc(10e-6)));
        c.addDevice(make_rcp<CurrentSource>("IB2", gnd, nb2,
                                            SourceWaveform::dc(10e-6)));
        c.addDevice(make_rcp<Bjt>("Q1", nc1, nb1, gnd, m));
        c.addDevice(make_rcp<Bjt>("Q2", nc2, nb2, gnd, m));
        DCAnalysis dc(c);
        DCResult r = dc.run(false, false);
        checkTrue("shared-card DC converges", r.success);
        if (r.success) {
            Bjt* q1 = static_cast<Bjt*>(c.getDevice("Q1"));
            Bjt* q2 = static_cast<Bjt*>(c.getDevice("Q2"));
            checkTrue("same card object shared",
                      q1->model().get() == q2->model().get());
            check("symmetric collector voltages",
                  nc1->getVoltage(), nc2->getVoltage(), 1e-9);
        }
    }

    // Netlist end-to-end (NPN + PNP in one file).
    {
        Circuit c("net");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL QN NPN (IS=1e-16 BF=100)\n"
            ".MODEL QP PNP (IS=1e-16 BF=50)\n"
            "VCC 1 0 DC 10\n"
            "RC 1 2 1k\n"
            "IB 0 3 DC 10u\n"
            "Q1 2 3 0 QN\n"
            "Q2 2 3 0 QP AREA=2\n"
            ".OP\n.END\n",
            c);
        checkTrue("netlist with NPN+PNP parses", ok);
        if (!ok) {
            std::printf("  [diag] parse error: %s\n",
                        p.getErrorMessage().c_str());
        }
        if (ok) {
            DCAnalysis dc(c);
            DCResult r = dc.run(false, false);
            checkTrue("netlist DC converges", r.success);
            if (r.success) {
                Bjt* q1 = static_cast<Bjt*>(c.getDevice("Q1"));
                Bjt* q2 = static_cast<Bjt*>(c.getDevice("Q2"));
                checkTrue("Q1 is NPN", !q1->model()->isPnp());
                checkTrue("Q2 is PNP", q2->model()->isPnp());
                check("Q2 AREA parsed", q2->instanceParams().AREA, 2.0, 0.0);
            }
        }
    }
    {
        Circuit c("e1");
        NetlistParser p;
        const bool ok = p.parseString("Q1 1 2 3 MISSING\n.END\n", c);
        checkTrue("unknown model rejected", !ok);
    }
    {
        Circuit c("e2");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-14)\nQ1 1 2 3 DM\n.END\n", c);
        checkTrue("diode card rejected for Q1", !ok);
    }
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    std::setvbuf(stdout, nullptr, _IONBF, 0);   // crash-safe diagnostics
    std::printf("=== BJT (Ebers-Moll) device verification ===\n");
    test_registry();
    test_aliases();
    test_active_region();
    test_early();
    test_saturation_and_cutoff();
    test_pnp_mirror();
    test_area();
    test_temperature();
    test_charge_model();
    test_transient_switching();
    test_sharing_and_parser();

    std::printf("\n=== %s (%d failure%s) ===\n",
                g_failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    mpi.finalize();
    return g_failures == 0 ? 0 : 1;
}
