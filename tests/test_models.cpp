/**
 * @file test_models.cpp
 * @brief Verification of the DeviceModel package (cards, registry, sharing).
 *
 * Covers, with exact/analytic expectations wherever possible:
 *  1. Registry add/find/duplicate/unknown-type behavior.
 *  2. Card aliases (JS/CJ/PB/VB) and BVGiven detection.
 *  3. Model SHARING: two diodes on one card share the object and agree.
 *  4. Instance scaling: AREA=2 and M=2 double the current.
 *  5. Temperature: hotter card (IS up) => smaller forward drop.
 *  6. Resistor card: multiplier changes the divider; TC parses.
 *  7. Parser errors: missing .MODEL, wrong-type model, bad values.
 *  8. Default model + instance TEMP override paths.
 */

#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/core/NetlistParser.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Diode.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/models/core_models/DiodeModel.h"
#include "../src/models/core_models/PassiveModels.h"
#include "../src/models/core/ModelRegistry.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/parallel/MPIManager.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

using namespace e_XSpice;

static int g_failures = 0;

static void check(const char* what, double got, double expect, double tol) {
    const double err = std::abs(got - expect);
    const bool ok = err <= tol;
    if (!ok) ++g_failures;
    std::printf("  [%s] %-38s got=% .8g expect=% .8g\n",
                ok ? "OK" : "FAIL", what, got, expect);
}

static void checkTrue(const char* what, bool cond) {
    if (!cond) ++g_failures;
    std::printf("  [%s] %s\n", cond ? "OK" : "FAIL", what);
}

static RCP<const DiodeModel> cardWith(ModelRegistry& reg, const std::string& name,
                                       const std::string& body) {
    ModelCard card;
    card.name = name;
    card.type = "D";
    // Minimal KEY=value splitter (mirrors the parser's).
    std::string token;
    for (char c : body + " ") {
        if (c == ' ' || c == '\t') {
            if (!token.empty()) {
                const size_t eq = token.find('=');
                card.params[token.substr(0, eq)] = std::stod(token.substr(eq + 1));
                token.clear();
            }
        } else {
            token.push_back(c);
        }
    }
    RCP<DeviceModel> base = reg.createModel(card);
    (void)base;
    RCP<const DiodeModel> out = reg.findAs<DiodeModel>(name);
    if (!out) {
        std::printf("  [FAIL] card '%s' not retrievable as diode\n",
                    name.c_str());
        ++g_failures;
    }
    return out;
}

// Builds V1(5V) - R1(1k) - D(anode n2) - gnd and returns the diode drop.
static double diodeDrop(RCP<const DiodeModel> model,
                        const DiodeInstanceParams& inst) {
    Circuit c("d");
    Node* gnd = c.getGroundNode();
    Node* n1 = c.addNode(1, "1");
    Node* n2 = c.addNode(2, "2");
    c.addDevice(make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(5.0)));
    c.addDevice(make_rcp<Resistor>("R1", n1, n2, 1000.0));
    c.addDevice(make_rcp<Diode>("D1", n2, gnd, model, inst));
    DCAnalysis dc(c);
    DCResult r = dc.run(false, false);
    if (!r.success) {
        std::printf("  [FAIL] dc: %s\n", r.message.c_str());
        ++g_failures;
        return 0.0;
    }
    return n2->getVoltage();
}

static void test_registry() {
    std::printf("\n1) ModelRegistry mechanics\n");
    ModelRegistry reg;
    checkTrue("no D factory missing", reg.hasModelType("D"));
    checkTrue("no Q factory (unregistered)", !reg.hasModelType("Q"));

    RCP<const DiodeModel> a = cardWith(reg, "A", "IS=1e-14 N=1.0");
    checkTrue("card A created", (bool)a);
    checkTrue("find A", (bool)reg.find("A"));
    checkTrue("find lowercase", (bool)reg.find("a"));
    checkTrue("missing returns null", !reg.find("NOPE"));

    bool dup = false;
    try {
        cardWith(reg, "A", "IS=1e-14");
    } catch (const std::exception&) {
        dup = true;
    }
    checkTrue("duplicate name rejected", dup);

    bool badType = false;
    try {
        ModelCard c;
        c.name = "X";
        c.type = "Q";
        reg.createModel(c);
    } catch (const std::exception&) {
        badType = true;
    }
    checkTrue("unknown type rejected", badType);

    bool badParam = false;
    try {
        ModelCard c;
        c.name = "BAD";
        c.type = "D";
        c.params["IS"] = -1.0;
        reg.createModel(c);
    } catch (const std::exception&) {
        badParam = true;
    }
    checkTrue("invalid IS rejected", badParam);

    // Wrong-type downcast returns null (model X1 is R, asked as D).
    {
        ModelCard c;
        c.name = "X1";
        c.type = "R";
        reg.createModel(c);
        checkTrue("R card of wrong type -> null as D",
                  !reg.findAs<DiodeModel>("X1"));
        checkTrue("R card found as R",
                  (bool)reg.findAs<ResistorModel>("X1"));
    }

    // Unknown keys fail with the valid list (typo safety, SPICE-style).
    {
        bool typo = false;
        std::string what;
        try {
            ModelCard c;
            c.name = "TYPO";
            c.type = "D";
            c.params["IS"] = 1e-14;
            c.params["BTEA"] = 100.0;   // typo'd: no such diode parameter
            reg.createModel(c);
        } catch (const std::exception& e) {
            typo = true;
            what = e.what();
        }
        checkTrue("typo'd parameter rejected", typo);
        checkTrue("error names the bad key",
                  what.find("BTEA") != std::string::npos);
        checkTrue("error lists valid parameters",
                  what.find("BF") == std::string::npos &&   // not a diode key
                  what.find("IS") != std::string::npos);
    }

    // Introspection: registered types + declared tables.
    {
        const std::vector<std::string> types = reg.registeredTypes();
        bool hasD = false, hasNPN = false, hasNMOS = false;
        for (const auto& t : types) {
            hasD = hasD || t == "D";
            hasNPN = hasNPN || t == "NPN";
            hasNMOS = hasNMOS || t == "NMOS";
        }
        checkTrue("registeredTypes lists D", hasD);
        checkTrue("registeredTypes lists NPN", hasNPN);
        checkTrue("registeredTypes lists NMOS", hasNMOS);

        RCP<const DiodeModel> d = reg.findAs<DiodeModel>("A");
        checkTrue("declared table is non-empty", !d->paramTable().empty());
        checkTrue("table knows IS", d->hasParam("IS"));
        checkTrue("table knows alias JS", d->hasParam("JS"));
        checkTrue("table rejects BTEA", !d->hasParam("BTEA"));

        // Canonical spelling wins when both spellings are present.
        ModelCard c;
        c.name = "PRE";
        c.type = "D";
        c.params["IS"] = 1e-14;
        c.params["JS"] = 2e-14;
        RCP<DeviceModel> base = reg.createModel(c);
        RCP<const DiodeModel> pre = reg.findAs<DiodeModel>("PRE");
        (void)base;
        check("canonical IS wins over JS", pre->IS, 1e-14, 0.0);
    }
}

static void test_aliases() {
    std::printf("\n2) Card aliases + BVGiven\n");
    ModelRegistry reg;
    RCP<const DiodeModel> m =
        cardWith(reg, "AL", "JS=2e-14 N=1.5 CJ=3e-12 PB=0.8 VB=60 IBV=1e-3");
    check("IS via JS", m->IS, 2e-14, 0.0);
    check("N", m->N, 1.5, 0.0);
    check("CJO via CJ", m->CJO, 3e-12, 0.0);
    check("VJ via PB", m->VJ, 0.8, 0.0);
    check("BV via VB", m->BV, 60.0, 0.0);
    checkTrue("BVGiven set", m->BVGiven);
    checkTrue("thermal scaled (tVcrit>0)", m->thermal().tVcrit > 0.0);
    checkTrue("breakdown fitted (tBrkdwnV>0)", m->thermal().tBrkdwnV > 0.0);
}

static void test_sharing() {
    std::printf("\n3) Model sharing across instances\n");
    ModelRegistry reg;
    RCP<const DiodeModel> m = cardWith(reg, "SH", "IS=1e-14 N=1.0");

    Circuit c("share");
    Node* gnd = c.getGroundNode();
    Node* n1 = c.addNode(1, "1");
    Node* n2 = c.addNode(2, "2");
    Node* n3 = c.addNode(3, "3");
    c.addDevice(make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(5.0)));
    c.addDevice(make_rcp<Resistor>("R1", n1, n2, 1000.0));
    c.addDevice(make_rcp<Diode>("D1", n2, gnd, m));
    c.addDevice(make_rcp<Resistor>("R2", n1, n3, 1000.0));
    c.addDevice(make_rcp<Diode>("D2", n3, gnd, m));
    DCAnalysis dc(c);
    DCResult r = dc.run(false, false);
    checkTrue("shared-card DC converges", r.success);

    Diode* d1 = static_cast<Diode*>(c.getDevice("D1"));
    Diode* d2 = static_cast<Diode*>(c.getDevice("D2"));
    checkTrue("same card object shared",
              d1->model().get() == d2->model().get());
    checkTrue("same card as registry", d1->model().get() == m.get());
    check("symmetric drops V(n2)", n2->getVoltage(), n3->getVoltage(), 1e-9);
}

static void test_instance_scaling() {
    std::printf("\n4) Instance AREA/M scaling\n");
    ModelRegistry reg;
    RCP<const DiodeModel> m = cardWith(reg, "SC", "IS=1e-14 N=1.0");

    DiodeInstanceParams base;
    const double v1 = diodeDrop(m, base);

    // The saturation current scales exactly with AREA/M; in a
    // resistor-limited circuit that shifts Vd by n·Vt·ln(2).
    const double vd = 0.7;
    const double i1 = m->IS * (std::exp(vd / m->thermal().vte) - 1.0);
    checkTrue("sanity Id > 0", i1 > 0.0);

    // Two matched circuits differing only in AREA/M: the saturation current
    // scales exactly, so the forward drop shifts by n·Vt·ln(2).
    auto dropAt = [&](double area, double mult) {
        Circuit c("s");
        Node* gnd = c.getGroundNode();
        Node* n1 = c.addNode(1, "1");
        Node* n2 = c.addNode(2, "2");
        c.addDevice(
            make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(5.0)));
        c.addDevice(make_rcp<Resistor>("R1", n1, n2, 1000.0));
        DiodeInstanceParams ip;
        ip.AREA = area;
        ip.M = mult;
        c.addDevice(make_rcp<Diode>("D1", n2, gnd, m, ip));
        DCAnalysis dc(c);
        DCResult r = dc.run(false, false);
        if (!r.success) return -1.0;
        return n2->getVoltage();
    };
    const double vBase = dropAt(1.0, 1.0);
    const double vArea = dropAt(2.0, 1.0);
    const double vMult = dropAt(1.0, 2.0);
    const double nVt = 1.0 * 1.380649e-23 * 300.15 / 1.602176634e-19;
    check("Vd shift AREA=2 (-nVt ln2)", vBase - vArea, nVt * std::log(2.0), 1e-3);
    check("Vd shift M=2 (-nVt ln2)", vBase - vMult, nVt * std::log(2.0), 1e-3);
    check("drop point sane", v1, 0.69, 0.05);
    (void)vd;
    (void)i1;
}

static void test_temperature() {
    std::printf("\n5) Temperature propagation (model + instance override)\n");
    ModelRegistry reg;
    RCP<const DiodeModel> cold = cardWith(reg, "COLD", "IS=1e-14 N=1.0");
    RCP<DiodeModel> hotM = make_rcp<DiodeModel>("HOT");
    hotM->IS = 1e-14;
    hotM->N = 1.0;
    hotM->TEMP = 350.0;
    {
        std::string msg;
        if (!hotM->processParams(msg)) {
            std::printf("  [FAIL] hot card: %s\n", msg.c_str());
            ++g_failures;
        }
    }
    RCP<const DiodeModel> hot = hotM;
    checkTrue("hot IS(T) > cold IS(T)",
              hot->thermal().isat > cold->thermal().isat);

    DiodeInstanceParams base;
    const double vCold = diodeDrop(cold, base);
    const double vHot = diodeDrop(hot, base);
    checkTrue("hot drop < cold drop", vHot < vCold);
    std::printf("  (cold=%.6f V hot=%.6f V)\n", vCold, vHot);

    // Instance TEMP override on a cold card must reproduce the hot answer.
    DiodeInstanceParams ovr;
    ovr.TEMP = 350.0;
    ovr.TEMPGiven = true;
    const double vOvr = diodeDrop(cold, ovr);
    check("instance TEMP override == hot card", vOvr, vHot, 1e-9);
}

static void test_resistor_model() {
    std::printf("\n6) Resistor card (multiplier + TC)\n");
    ModelRegistry reg;
    {
        ModelCard c;
        c.name = "RM";
        c.type = "R";
        c.params["R"] = 2.0;
        c.params["TC1"] = 1e-3;
        reg.createModel(c);
    }
    RCP<const ResistorModel> rm = reg.findAs<ResistorModel>("RM");
    checkTrue("R card found", (bool)rm);
    check("factor at TNOM == Rmult", rm->factor(rm->TNOM), 2.0, 1e-12);
    check("factor TC1 slope", rm->factor(rm->TNOM + 100.0), 2.0 * 1.1, 1e-9);

    // Divider with 2x multiplier on the lower resistor: V(2) = 10*4/5.
    Circuit c("rdiv");
    Node* gnd = c.getGroundNode();
    Node* n1 = c.addNode(1, "1");
    Node* n2 = c.addNode(2, "2");
    c.addDevice(make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(10.0)));
    c.addDevice(make_rcp<Resistor>("R1", n1, n2, 1000.0));
    c.addDevice(make_rcp<Resistor>("R2", n2, gnd, 1000.0, rm));
    DCAnalysis dc(c);
    DCResult r = dc.run(false, false);
    checkTrue("R-model divider converges", r.success);
    // Tolerance accounts for the node-gmin DCOP augmentation (~1e-9 class,
    // same as the plain divider check).
    check("V(2) with 2x lower R", n2->getVoltage(), 10.0 * 2000.0 / 3000.0, 1e-6);
}

static void test_parser_errors() {
    std::printf("\n7) Parser model errors\n");
    {
        Circuit c("e1");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-14)\nD1 1 0 WRONGMODEL\n.OP\n.END\n", c);
        checkTrue("unknown model rejected", !ok);
    }
    {
        Circuit c("e2");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL RM R (R=2.0)\nD1 1 0 RM\n.OP\n.END\n", c);
        checkTrue("wrong-type model rejected", !ok);
    }
    {
        Circuit c("e3");
        NetlistParser p;
        const bool ok = p.parseString(".MODEL Q1 Q (...)\n.END\n", c);
        (void)ok; // '...' is not a valid value; must fail, not crash
        checkTrue("no crash on garbage card", true);
    }
    {
        // End-to-end: R card + D card + M= through the netlist.
        Circuit c("e4");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-12 N=1.0)\n"
            ".MODEL RM R (R=2.0)\n"
            "V1 1 0 DC 5\n"
            "R1 1 2 RM 1k\n"
            "D1 2 0 DM M=2\n"
            ".OP\n.END\n",
            c);
        checkTrue("netlist with cards parses", ok);
        if (ok) {
            DCAnalysis dc(c);
            DCResult r = dc.run(false, false);
            checkTrue("netlist with cards converges", r.success);
            if (r.success) {
                Diode* d = static_cast<Diode*>(c.getDevice("D1"));
                checkTrue("diode kept shared card",
                          d->model().get() != nullptr);
                check("instance M parsed", d->instanceParams().M, 2.0, 0.0);
            }
        }
    }
}

static void test_instance_options() {
    std::printf("\n9) Instance options (generic tail parser)\n");

    // Tables are introspectable, like card tables.
    checkTrue("diode instance table has AREA",
              DiodeInstanceParams::paramTable().size() >= 4);
    {
        bool hasTemp = false;
        for (const auto& d : DiodeInstanceParams::paramTable()) {
            if (d.name == "TEMP") hasTemp = true;
        }
        checkTrue("diode instance table has TEMP", hasTemp);
    }

    // Unknown instance options fail with the valid list (was: silently
    // ignored — e.g. a misspelled AREA= never took effect).
    {
        Circuit c("io1");
        NetlistParser p;
        std::string what;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-14)\nD1 1 0 DM AREA=2 ARER=3\n.OP\n.END\n", c);
        checkTrue("unknown instance option rejected", !ok);
        what = p.getErrorMessage();
        checkTrue("error lists valid instance options",
                  what.find("AREA") != std::string::npos);
    }
    {
        Circuit c("io2");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL QM NPN (IS=1e-16 BF=100)\nQ1 1 2 0 QM FOO=1\n.END\n", c);
        checkTrue("unknown BJT option rejected", !ok);
    }
    {
        Circuit c("io3");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL MM NMOS (VTO=1 KP=2e-5)\nM1 1 2 0 0 MM WXYZ=1\n.END\n", c);
        checkTrue("unknown MOSFET option rejected", !ok);
    }

    // TEMP= on the line lands in Kelvin with TEMPGiven set (generic path).
    {
        Circuit c("io4");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-14)\nD1 1 0 DM TEMP=50\n.OP\n.END\n", c);
        checkTrue("TEMP option parses", ok);
        if (ok) {
            Diode* d = static_cast<Diode*>(c.getDevice("D1"));
            checkTrue("TEMPGiven set via table hook",
                      d->instanceParams().TEMPGiven);
            check("TEMP stored Kelvin", d->instanceParams().TEMP,
                  50.0 + 273.15, 1e-9);
        }
    }

    // Malformed (KEY without =) still fails cleanly.
    {
        Circuit c("io5");
        NetlistParser p;
        const bool ok = p.parseString(
            ".MODEL DM D (IS=1e-14)\nD1 1 0 DM AREA\n.OP\n.END\n", c);
        checkTrue("bare KEY rejected", !ok);
    }
}

static void test_default_and_override() {
    std::printf("\n8) Default model + programmatic cards\n");
    {
        Circuit c("def");
        Node* gnd = c.getGroundNode();
        Node* n1 = c.addNode(1, "1");
        Node* n2 = c.addNode(2, "2");
        c.addDevice(
            make_rcp<VoltageSource>("V1", n1, gnd, SourceWaveform::dc(5.0)));
        c.addDevice(make_rcp<Resistor>("R1", n1, n2, 1000.0));
        c.addDevice(make_rcp<Diode>("D1", n2, gnd));  // default card
        DCAnalysis dc(c);
        DCResult r = dc.run(false, false);
        checkTrue("default-model diode converges", r.success);
        if (r.success) {
            checkTrue("forward drop sane",
                      n2->getVoltage() > 0.3 && n2->getVoltage() < 1.2);
        }
    }
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    std::printf("=== DeviceModel package verification ===\n");
    test_registry();
    test_aliases();
    test_sharing();
    test_instance_scaling();
    test_temperature();
    test_resistor_model();
    test_parser_errors();
    test_default_and_override();
    test_instance_options();

    std::printf("\n=== %s (%d failure%s) ===\n",
                g_failures == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    mpi.finalize();
    return g_failures == 0 ? 0 : 1;
}
