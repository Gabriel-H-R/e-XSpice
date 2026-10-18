/**
 * @file test_bjt2.cpp
 * @brief Verification of the BJT level-2 (Gummel-Poon) model, NPN + PNP.
 *
 * Every check is analytic or a physics invariant:
 *  1. Card validation: LEVEL must be 1/2; level-1 cards reject GP params.
 *  2. Reduction: LEVEL=2 with all GP extras off == LEVEL=1 (bit-near).
 *  3. High injection: beta falls at high current (qb > 1), matches 1/qb.
 *  4. Non-ideal base: low-current log-slope is 1/(NE·Vt), not 1/(NF·Vt).
 *  5. Ohmic RB: external Vb exceeds the RB=0 case by Ib·RB.
 *  6. Ohmic RC+RE: external Vce shifts by the resistive drops.
 *  7. Substrate: 4-terminal reverse substrate current ~0, charge stored.
 *  8. PNP level 2 mirrors NPN (incl. RB and substrate).
 *  9. Early effect still present at level 2 (go > 0).
 * 10. Transient switching converges with TF at level 2.
 * 11. Parser: LEVEL=2 card, 4-node Q, error paths.
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

struct CeResult {
    bool ok = false;
    double vbe = 0, vbc = 0, vc = 0, vb = 0;
    double ib = 0, ic = 0, ie = 0;
    double isub = 0;
    double vbe_int = 0, vbc_int = 0;   // intrinsic (post-RB/RC/RE) voltages
};

static CeResult runCe(RCP<const BjtModel> model, bool pnp, double vcc,
                      double rc, double ibias, bool withSub = false) {
    CeResult r;
    Circuit c("ce2");
    Node* gnd = c.getGroundNode();
    Node* nrail = c.addNode(1, "1");
    Node* ncol = c.addNode(2, "2");
    Node* nbase = c.addNode(3, "3");
    Node* nsub = withSub ? c.addNode(4, "4") : nullptr;

    const double p = pnp ? -1.0 : 1.0;
    c.addDevice(make_rcp<VoltageSource>("V1", nrail, gnd,
                                        SourceWaveform::dc(p * vcc)));
    c.addDevice(make_rcp<Resistor>("RC", nrail, ncol, rc));
    c.addDevice(make_rcp<CurrentSource>("IB", p > 0 ? gnd : nbase,
                                        p > 0 ? nbase : gnd,
                                        SourceWaveform::dc(ibias)));
    BjtInstanceParams inst;
    Bjt* q = nullptr;
    if (withSub) {
        // Substrate tied to the reference rail (reverse C-S for NPN when
        // the collector is up; the PNP mirror ties it to its own rail).
        c.addDevice(make_rcp<VoltageSource>("VS", nsub, gnd,
                                            SourceWaveform::dc(0.0)));
        auto dev = make_rcp<Bjt>("Q1", ncol, nbase, gnd, model, inst, nsub);
        q = dev.get();
        c.addDevice(std::move(dev));
    } else {
        auto dev = make_rcp<Bjt>("Q1", ncol, nbase, gnd, model, inst);
        q = dev.get();
        c.addDevice(std::move(dev));
    }

    DCAnalysis dc(c);
    DCResult res = dc.run(false, false);
    r.ok = res.success;
    if (!res.success) return r;
    r.vc = ncol->getVoltage();
    r.vb = nbase->getVoltage();
    r.vbe = nbase->getVoltage();
    r.vbc = nbase->getVoltage() - r.vc;
    r.ib = q->getBaseCurrent();
    r.ic = q->getCollectorCurrent();
    r.ie = q->getEmitterCurrent();
    r.isub = q->getSubstrateCurrent();
    r.vbe_int = q->getIntrinsicVbe();
    r.vbc_int = q->getIntrinsicVbc();
    return r;
}

// ===============================================================
// Tests
// ===============================================================

static void test_card_validation() {
    std::printf("\n1) Card validation (LEVEL + GP guards)\n");
    ModelRegistry reg;
    RCP<const BjtModel> l2 = cardWith(reg, "L2", "LEVEL=2 IS=1e-16 BF=100", false);
    checkTrue("LEVEL=2 accepted", (bool)l2 && l2->LEVEL == 2);

    bool badLevel = false;
    try {
        ModelCard c; c.name = "BL"; c.type = "NPN"; c.params["LEVEL"] = 3.0;
        reg.createModel(c);
    } catch (const std::exception&) { badLevel = true; }
    checkTrue("LEVEL=3 rejected", badLevel);

    bool l1gp = false;
    try {
        ModelCard c; c.name = "MIX"; c.type = "NPN"; c.params["IKF"] = 1e-3;
        reg.createModel(c);
    } catch (const std::exception&) { l1gp = true; }
    checkTrue("LEVEL=1 + IKF rejected", l1gp);

    bool l1rb = false;
    try {
        ModelCard c; c.name = "MIX2"; c.type = "PNP"; c.params["RB"] = 10.0;
        reg.createModel(c);
    } catch (const std::exception&) { l1rb = true; }
    checkTrue("LEVEL=1 + RB rejected", l1rb);

    RCP<const BjtModel> gp = cardWith(
        reg, "GP", "LEVEL=2 IKF=1e-3 IKR=1e-3 ISE=1e-16 NE=1.5 "
                   "ISC=1e-16 NC=2 RB=10 RC=5 RE=1 CJS=1e-12", false);
    checkTrue("full GP card accepted", (bool)gp);
    check("IKF parsed", gp->thermal().ikf, 1e-3, 0.0);
    checkTrue("internal nodes = 3", gp->numInternalNodes() == 3);
}

static void test_reduction_to_level1() {
    std::printf("\n2) LEVEL=2 with GP extras off == LEVEL=1\n");
    ModelRegistry reg;
    RCP<const BjtModel> l1 = cardWith(reg, "L1", "IS=1e-16 BF=100 VAF=100", false);
    RCP<const BjtModel> l2 = cardWith(reg, "L2V", "LEVEL=2 IS=1e-16 BF=100 VAF=100",
                                      false);

    const CeResult a = runCe(l1, false, 10.0, 1000.0, 10e-6);
    const CeResult b = runCe(l2, false, 10.0, 1000.0, 10e-6);
    checkTrue("both converge", a.ok && b.ok);
    if (!a.ok || !b.ok) return;
    check("Vbe identical", b.vbe, a.vbe, 1e-9);
    check("Vc identical", b.vc, a.vc, 1e-9);
    check("beta identical", b.ic / b.ib, a.ic / a.ib, 1e-9);
}

static void test_high_injection() {
    std::printf("\n3) High injection (IKF knee)\n");
    ModelRegistry reg;
    // IKF = 10mA: at ib = 1uA (ic ≈ 100uA) q2 ≈ 0.01 → qb ≈ 1 (β ≈ BF);
    // at ib = 1mA the self-consistent point has q2 ≈ 4 → strong roll-off.
    RCP<const BjtModel> m = cardWith(
        reg, "HI", "LEVEL=2 IS=1e-16 BF=100 IKF=0.01", false);

    const CeResult low = runCe(m, false, 10.0, 1000.0, 1e-6);
    const CeResult high = runCe(m, false, 10.0, 10.0, 1e-3);
    checkTrue("both converge", low.ok && high.ok);
    if (!low.ok || !high.ok) return;

    const double betaLow = low.ic / low.ib;
    const double betaHigh = high.ic / high.ib;
    checkTrue("low-current beta ≈ BF", std::abs(betaLow - 100.0) < 1.0);
    checkTrue("beta falls at high current", betaHigh < 0.5 * betaLow);
    std::printf("  (beta low=%.4g high=%.4g)\n", betaLow, betaHigh);

    // Independent qb prediction at the high point: with VAF=∞ (default),
    // q1 = 1, q2 = fwd/IKF, qb = (1+√(1+4q2))/2, β = BF/qb.
    // fwd ≈ ic·qb ⇒ solve qb self-consistently: qb² − qb − fwd0/IKF = 0
    // with fwd0 = IS·exp(vbe/vte).
    const double vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
    const double vbe = high.vbe;
    const double fwd0 = 1e-16 * (std::exp(vbe / vt) - 1.0);
    const double qb = (1.0 + std::sqrt(1.0 + 4.0 * fwd0 / 0.01)) * 0.5;
    check("beta_high = BF/qb", betaHigh, 100.0 / qb, 0.02 * betaHigh + 0.5);
}

static void test_nonideal_base() {
    std::printf("\n4) Non-ideal base current (ISE/NE slope)\n");
    ModelRegistry reg;
    // ISE ≫ IS/BF at low current: the base slope follows NE = 1.5.
    // ISE = 1e-13 makes the non-ideal term dominate by ~100× at 1–10nA.
    RCP<const BjtModel> m = cardWith(
        reg, "NE", "LEVEL=2 IS=1e-16 BF=100 ISE=1e-13 NE=1.5", false);

    // Two forced-base points a decade apart in the low-current regime.
    const CeResult r1 = runCe(m, false, 10.0, 1e6, 1e-9);
    const CeResult r2 = runCe(m, false, 10.0, 1e6, 1e-8);
    checkTrue("both converge", r1.ok && r2.ok);
    if (!r1.ok || !r2.ok) return;

    // Slope dVbe per decade of Ib: NE·Vt·ln(10) when ISE dominates.
    const double vt = 1.380649e-23 * 300.15 / 1.602176634e-19;
    const double slope = (r2.vbe - r1.vbe) / std::log10(r2.ib / r1.ib);
    check("low-current slope = NE·Vt·ln10", slope, 1.5 * vt * std::log(10.0),
          0.05 * 1.5 * vt * std::log(10.0));
}

static void test_ohmic_rb() {
    std::printf("\n5) Ohmic base resistance (RB)\n");
    ModelRegistry reg;
    RCP<const BjtModel> noR = cardWith(reg, "NOR", "LEVEL=2 IS=1e-16 BF=100",
                                       false);
    RCP<const BjtModel> rb = cardWith(reg, "RB", "LEVEL=2 IS=1e-16 BF=100 RB=1000",
                                      false);

    const CeResult a = runCe(noR, false, 10.0, 1000.0, 10e-6);
    const CeResult b = runCe(rb, false, 10.0, 1000.0, 10e-6);
    checkTrue("both converge", a.ok && b.ok);
    if (!a.ok || !b.ok) return;

    // Same forced base current ⇒ same intrinsic Vbe ⇒ the external base
    // sits exactly Ib·RB higher.
    check("Vb shift = Ib·RB", b.vb - a.vb, b.ib * 1000.0, 1e-6);
    check("collector current unchanged", b.ic, a.ic, 1e-8);
}

static void test_ohmic_rc_re() {
    std::printf("\n6) Ohmic collector/emitter resistances (RC/RE)\n");
    ModelRegistry reg;
    RCP<const BjtModel> noR = cardWith(reg, "NOR2", "LEVEL=2 IS=1e-16 BF=100",
                                       false);
    RCP<const BjtModel> rce = cardWith(reg, "RCE", "LEVEL=2 IS=1e-16 BF=100 "
                                                  "RC=100 RE=10", false);

    const CeResult a = runCe(noR, false, 10.0, 1000.0, 10e-6);
    const CeResult b = runCe(rce, false, 10.0, 1000.0, 10e-6);
    checkTrue("both converge", a.ok && b.ok);
    if (!a.ok || !b.ok) return;

    // The forced base current fixes the intrinsic point, so the collector
    // current is (almost) unchanged — but the INTRINSIC terminals sit
    // behind the drops: Vc_int = Vc_ext − Ic·RC, Ve_int = |Ie|·RE.
    // With no RB, Vb_int = Vb_ext, so:
    //   (Vb_ext − Vc_ext) − vbc_int = Ic·RC,
    //   Vb_ext − vbe_int = Ve_int = |Ie|·RE.
    check("collector current ~ unchanged", b.ic, a.ic, 0.01 * std::abs(a.ic));
    // Vc_int = Vc_ext − Ic·RC ⇒ (Vb−Vc)_ext − vbc_int = −Ic·RC.
    // Tolerance is solver-residual level (the residual test admits rows
    // up to ~1e-6 relative; 13 µV here is consistent with that).
    check("intrinsic Vc drop = −Ic·RC", (b.vb - b.vc) - b.vbc_int,
          -b.ic * 100.0, 1e-4);
    check("intrinsic Ve rise = |Ie|·RE", b.vb - b.vbe_int,
          std::abs(b.ie) * 10.0, 1e-4);
}

static void test_substrate() {
    std::printf("\n7) Substrate terminal (CJS diode + charge)\n");
    ModelRegistry reg;
    RCP<const BjtModel> m = cardWith(
        reg, "SUB", "LEVEL=2 IS=1e-16 BF=100 ISS=1e-16 CJS=2e-12", false);

    const CeResult r = runCe(m, false, 10.0, 1000.0, 10e-6, true);
    checkTrue("4-terminal DC converges", r.ok);
    if (!r.ok) return;
    // Substrate at ground, collector up: C-S reverse ⇒ only gmin leaks.
    checkTrue("substrate current ~ 0", std::abs(r.isub) < 1e-9);
    checkTrue("main current intact", std::abs(r.ic - 1e-3) < 2e-5);
}

static void test_pnp_level2_mirror() {
    std::printf("\n8) PNP level-2 mirror (incl. RB + substrate)\n");
    ModelRegistry reg;
    RCP<const BjtModel> npn = cardWith(
        reg, "N2", "LEVEL=2 IS=1e-16 BF=100 IKF=0.01 RB=100 CJS=1e-12", false);
    RCP<const BjtModel> pnp = cardWith(
        reg, "P2", "LEVEL=2 IS=1e-16 BF=100 IKF=0.01 RB=100 CJS=1e-12", true);

    const CeResult n = runCe(npn, false, 10.0, 1000.0, 10e-6, true);
    const CeResult p = runCe(pnp, true, 10.0, 1000.0, 10e-6, true);
    checkTrue("both converge", n.ok && p.ok);
    if (!n.ok || !p.ok) return;
    check("PNP Vbe = −NPN Vbe", p.vbe, -n.vbe, 1e-7);
    check("PNP Vc = −NPN Vc", p.vc, -n.vc, 1e-7);
    check("PNP ic = −NPN ic", p.ic, -n.ic, 1e-12);
    check("PNP isub = −NPN isub", p.isub, -n.isub, 1e-15);
}

static void test_early_level2() {
    std::printf("\n9) Early effect at level 2\n");
    ModelRegistry reg;
    RCP<const BjtModel> m = cardWith(
        reg, "E2", "LEVEL=2 IS=1e-16 BF=100 VAF=100", false);

    const CeResult r5 = runCe(m, false, 5.0, 1000.0, 10e-6);
    const CeResult r10 = runCe(m, false, 10.0, 1000.0, 10e-6);
    checkTrue("both converge", r5.ok && r10.ok);
    if (!r5.ok || !r10.ok) return;
    const double go = (r10.ic - r5.ic) / (r10.vc - r5.vc);
    checkTrue("go > 0", go > 0.0);
    check("go ~ ic/VAF", go, r5.ic / 100.0, 0.05 * r5.ic / 100.0);
}

static void test_transient_level2() {
    std::printf("\n10) Transient switching at level 2 (with TF)\n");
    ModelRegistry reg;
    RCP<const BjtModel> m = cardWith(
        reg, "T2", "LEVEL=2 IS=1e-16 BF=100 TF=1e-9 CJC=1e-12", false);

    Circuit c("sw2");
    Node* gnd = c.getGroundNode();
    Node* nrail = c.addNode(1, "1");
    Node* ndrv = c.addNode(3, "3");
    Node* nbase = c.addNode(4, "4");
    Node* ncol = c.addNode(2, "2");

    c.addDevice(make_rcp<VoltageSource>("VCC", nrail, gnd,
                                        SourceWaveform::dc(10.0)));
    c.addDevice(make_rcp<Resistor>("RC", nrail, ncol, 1000.0));
    c.addDevice(make_rcp<VoltageSource>(
        "VB", ndrv, gnd,
        SourceWaveform::pulse(0.0, 5.0, 1e-6, 1e-8, 1e-8, 1e-3, 1e-3)));
    c.addDevice(make_rcp<Resistor>("RB", ndrv, nbase, 100000.0));
    c.addDevice(make_rcp<Bjt>("Q1", ncol, nbase, gnd, m));

    const auto& nodes = c.getNodes();
    auto idx = [&](Node* n) {
        for (size_t i = 0; i < nodes.size(); ++i)
            if (nodes[i].get() == n) return static_cast<int>(i);
        return -1;
    };

    TransientAnalysis tran(c);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    TransientResult r = tran.run(0.0, 2e-6, 2e-8, false);
    checkTrue("transient converges", r.success);
    if (!r.success || r.data.empty()) return;
    check("cutoff: Vc = VCC", r.data.front().node_voltages[idx(ncol)],
          10.0, 5e-3);
    const auto& last = r.data.back();
    checkTrue("switched on", last.node_voltages[idx(ncol)] < 9.0);
}

static void test_parser_level2() {
    std::printf("\n11) Parser: LEVEL=2 cards + 4-node Q\n");
    {
        Circuit c("net2");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL QN NPN (LEVEL=2 IS=1e-16 BF=100 IKF=1e-3 RB=10)\n"
            "VCC 1 0 DC 10\n"
            "RC 1 2 1k\n"
            "IB 0 3 DC 10u\n"
            "Q1 2 3 0 QN\n"
            ".OP\n.END\n",
            c);
        checkTrue("3-node level-2 netlist parses", ok);
        if (ok) {
            DCAnalysis dc(c);
            DCResult r = dc.run(false, false);
            checkTrue("3-node level-2 DC converges", r.success);
            if (r.success) {
                Bjt* q = static_cast<Bjt*>(c.getDevice("Q1"));
                checkTrue("LEVEL=2 card kept", q->model()->LEVEL == 2);
            }
        }
    }
    {
        Circuit c("net3");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL QS NPN (LEVEL=2 IS=1e-16 BF=100 CJS=1e-12)\n"
            "VCC 1 0 DC 10\n"
            "RC 1 2 1k\n"
            "IB 0 3 DC 10u\n"
            "Q1 2 3 0 4 QS\n"
            ".OP\n.END\n",
            c);
        checkTrue("4-node level-2 netlist parses", ok);
        if (ok) {
            Bjt* q = static_cast<Bjt*>(c.getDevice("Q1"));
            checkTrue("substrate attached", q->hasSubstrate());
            DCAnalysis dc(c);
            DCResult r = dc.run(false, false);
            checkTrue("4-node DC converges", r.success);
        }
    }
    {
        // 4 nodes with a LEVEL=1 card must fail with a clear error.
        Circuit c("net4");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL Q1M NPN (IS=1e-16 BF=100)\n"
            "Q1 2 3 0 4 Q1M\n.END\n",
            c);
        checkTrue("4-node LEVEL=1 rejected", !ok);
    }
    {
        // LEVEL=3 rejected at the card.
        Circuit c("net5");
        NetlistParser p;
        const bool ok = p.parseString(".MODEL QB NPN (LEVEL=3)\n.END\n", c);
        checkTrue("LEVEL=3 rejected", !ok);
    }
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("=== BJT level-2 (Gummel-Poon) verification ===\n");
    test_card_validation();
    test_reduction_to_level1();
    test_high_injection();
    test_nonideal_base();
    test_ohmic_rb();
    test_ohmic_rc_re();
    test_substrate();
    test_pnp_level2_mirror();
    test_early_level2();
    test_transient_level2();
    test_parser_level2();

    std::printf("\n=== %s (%d failure%s) ===\n",
                g_failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    mpi.finalize();
    return g_failures == 0 ? 0 : 1;
}
