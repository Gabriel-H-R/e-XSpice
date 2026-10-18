/**
 * @file test_dc_sweep.cpp
 * @brief Verification of the DC sweep analysis (.DC).
 *
 * Every check is analytic or a structural invariant:
 *  1. Parser: `.DC [LIN|DEC|OCT] <var> <start> <stop> <step|points>` with an
 *     optional second (inner) specification, plus .OPTIONS (NOFAIL, RETRY,
 *     MAXITERS) and unknown-option warnings.
 *  2. Grid builder: LIN ascending/descending/single-point, DEC and OCT
 *     decade grids with the exact stop pinned.
 *  3. Sweep targets: V source, current source, resistor; precise errors for
 *     a missing device or a non-sweepable one.
 *  4. Resistor divider sweep (analytic Vout = Vin/2 at every point) with
 *     warm-start chaining; swept source restored afterwards.
 *  5. Nested sweep: two sources, both levels, analytic node voltages.
 *  6. Predictor equivalence: extrapolation must not change the solution.
 *  7. Failure ladder: conflicting parallel voltage sources have no DC
 *     solution -> abort by default, NOFAIL keeps sweeping, and the
 *     source-step refinement is exercised on the way.
 *  8. Streaming CSV: every point reaches the stream at full precision while
 *     the in-memory sample stays decimated (bounded RAM).
 *
 * Runs under plain ctest (size==1) and under mpiexec (same deterministic
 * sweep on every rank; distribution stays inside DaeSystem).
 */

#include "../src/core/Circuit.h"
#include "../src/core/Device.h"
#include "../src/core/Node.h"
#include "../src/core/NetlistParser.h"
#include "../src/devices/Diode.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/models/core_models/PassiveModels.h"
#include "../src/models/core_models/DiodeModel.h"
#include "../src/models/core/ModelRegistry.h"
#include "../src/analysis/DCAnalysis.h"
#include "../src/analysis/DCSweep.h"
#include "../src/parallel/MPIManager.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace e_XSpice;

static int g_failures = 0;

static void checkTrue(const char* what, bool cond) {
    if (!cond) ++g_failures;
    std::printf("  [%s] %s\n", cond ? "OK" : "FAIL", what);
}

static void checkClose(const char* what, double got, double expect,
                       double tol) {
    const bool ok = std::isfinite(got) && std::abs(got - expect) <= tol;
    if (!ok) ++g_failures;
    std::printf("  [%s] %-40s got=% .8g expect=% .8g\n",
                ok ? "OK" : "FAIL", what, got, expect);
}

// Divider: V1 -R1- out -R2- gnd (R1 == R2 by default).
static void buildDivider(Circuit& c, double vin, Node*& out,
                         double r1 = 1e3, double r2 = 1e3) {
    Node* gnd = c.getGroundNode();
    Node* top = c.addNode(1, "in");
    out = c.addNode(2, "out");
    c.addDevice(make_rcp<VoltageSource>("V1", top, gnd,
                                        SourceWaveform::dc(vin)));
    c.addDevice(make_rcp<Resistor>("R1", top, out, r1));
    c.addDevice(make_rcp<Resistor>("R2", out, gnd, r2));
}

static int nodeIndex(const Circuit& c, const Node* n) {
    const auto& nodes = c.getNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].get() == n) return static_cast<int>(i);
    }
    return -1;
}

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);
    const bool use_mpi = (mpi.getSize() > 1);
    const bool master = mpi.isMaster();
    bool ok = true;

    if (master) std::printf("DC sweep checks (size=%d)\n", mpi.getSize());

    // ---- 1) parser ----
    if (master) std::printf("\n1) .DC / .OPTIONS parsing\n");
    {
        {   // flat, LIN by default
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 in 0 DC 0\nR1 in out 1k\nR2 out 0 1k\n"
                ".DC VIN 0 10 1\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue("parses", parsed);
            checkTrue("want_dc_sweep set", p.want_dc_sweep);
            checkTrue("LIN default", p.dc_variation == "LIN");
            checkTrue("not nested", !p.dc_nested);
            checkTrue("source kept", p.dc_source == "VIN");
            checkClose("start", p.dc_start, 0.0, 0.0);
            checkClose("stop", p.dc_stop, 10.0, 0.0);
            checkClose("step", p.dc_step, 1.0, 0.0);
        }
        {   // DEC
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 in 0 DC 0\nR1 in out 1k\nR2 out 0 1k\n"
                ".DC DEC VIN 1 1000 2\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue("DEC parses", parsed);
            checkTrue("variation DEC", p.dc_variation == "DEC");
            checkClose("points-per-decade", p.dc_step, 2.0, 0.0);
        }
        {   // nested, second keyword explicit
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 a 0 DC 0\nV2 b 0 DC 0\nR1 a b 1k\n"
                ".DC LIN V1 0 2 1 DEC V2 1 100 3\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue("nested parses", parsed && p.dc_nested);
            checkTrue("outer source", p.dc_source == "V1");
            checkTrue("inner source", p.dc_source2 == "V2");
            checkTrue("inner variation DEC", p.dc_variation2 == "DEC");
            checkClose("inner start", p.dc_start2, 1.0, 0.0);
            checkClose("inner stop", p.dc_stop2, 100.0, 0.0);
            checkClose("inner points", p.dc_step2, 3.0, 0.0);
        }
        {   // nested without a keyword on the second spec (LIN default)
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 a 0 DC 0\nV2 b 0 DC 0\nR1 a b 1k\n"
                ".DC V1 0 2 1 V2 0 1 0.5\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue("nested (no keyword) parses", parsed && p.dc_nested);
            checkTrue("inner LIN default", p.dc_variation2 == "LIN");
            checkClose("inner step", p.dc_step2, 0.5, 0.0);
        }
        {   // three specifications -> rejected
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 a 0 DC 0\nV2 b 0 DC 0\nV3 c 0 DC 0\n"
                ".DC V1 0 2 1 V2 0 1 0.5 V3 0 1 0.5\n.END\n", c);
            checkTrue("three .DC specs rejected", !parsed);
        }
        {   // short card -> rejected
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 in 0 DC 0\nR1 in out 1k\nR2 out 0 1k\n"
                ".DC V1\n.END\n", c);
            checkTrue("short .DC rejected", !parsed);
        }
        {   // .OPTIONS
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 in 0 DC 0\nR1 in out 1k\nR2 out 0 1k\n"
                ".OPTIONS NOFAIL RETRY=5 MAXITERS=250 UNKNOWNIOP\n"
                ".DC V1 0 1 1\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue(".OPTIONS parses", parsed);
            checkTrue("NOFAIL", p.opt_nofail);
            checkTrue("RETRY=5", p.opt_max_retries == 5);
            checkTrue("MAXITERS=250", p.opt_max_iterations == 250);
            checkTrue("unknown option collected",
                      p.option_warnings.size() == 1
                      && p.option_warnings[0] == "UNKNOWNIOP");
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 2) grid builder ----
    if (master) std::printf("\n2) Sweep grids\n");
    {
        Circuit c("g");
        auto gridOf = [&c](DCSweepVariation v, double s, double e, double p) {
            DCSweep s2(c);
            s2.setVariation(v);
            s2.setRange(s, e, p);
            return s2.buildGrid();
        };
        {
            const auto g = gridOf(DCSweepVariation::LIN, 0.0, 10.0, 5.0);
            checkTrue("LIN point count", g.size() == 3);
            if (g.size() == 3) {
                checkClose("LIN grid", g[1], 5.0, 1e-12);
                checkClose("LIN stop pinned", g[2], 10.0, 1e-12);
            }
        }
        {
            const auto g = gridOf(DCSweepVariation::LIN, 10.0, 0.0, -5.0);
            checkTrue("LIN descending count", g.size() == 3);
            if (g.size() == 3) checkClose("LIN descending", g[1], 5.0, 1e-12);
        }
        {
            const auto g = gridOf(DCSweepVariation::LIN, 3.0, 3.0, 0.0);
            checkTrue("zero step = single point",
                      g.size() == 1 && g[0] == 3.0);
        }
        {
            const auto g = gridOf(DCSweepVariation::DEC, 1.0, 1000.0, 1.0);
            checkTrue("DEC 1/dec count", g.size() == 4);
            if (g.size() == 4) {
                checkClose("DEC decades", g[2], 100.0, 1e-9);
                checkClose("DEC stop pinned", g[3], 1000.0, 1e-9);
            }
        }
        {
            const auto g = gridOf(DCSweepVariation::OCT, 1.0, 8.0, 1.0);
            checkTrue("OCT 1/oct count", g.size() == 4);
            if (g.size() == 4) checkClose("OCT stop pinned", g[3], 8.0, 1e-9);
        }
        {
            const auto g = gridOf(DCSweepVariation::DEC, 100.0, 1.0, 2.0);
            checkTrue("degenerate DEC = single start point",
                      g.size() == 1 && g[0] == 100.0);
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 3) sweep targets ----
    if (master) std::printf("\n3) Sweepable targets\n");
    {
        Circuit c("t");
        Node* out = nullptr;
        buildDivider(c, 5.0, out);
        std::string err;

        auto t = DCSweep::makeTarget(c, "V1", err);
        checkTrue("voltage source target", t != nullptr);
        if (t) {
            t->set(3.5);
            checkClose("V source value", t->get(), 3.5, 0.0);
        }
        auto r = DCSweep::makeTarget(c, "R1", err);
        checkTrue("resistor target", r != nullptr);
        if (r) {
            const double original = r->get();
            checkClose("resistor initial value", original, 1e3, 1e-6);
            r->set(2e3);
            checkClose("resistor swept value", r->get(), 2e3, 1e-6);
            r->set(original);
        }
        checkTrue("missing device rejected",
                  DCSweep::makeTarget(c, "NOPE", err) == nullptr
                  && err.find("no device") != std::string::npos);
        // A diode is not sweepable: precise message, no crash.
        Node* gnd = c.getGroundNode();
        c.addDevice(make_rcp<Diode>("D1", out, gnd));
        err.clear();
        checkTrue("non-sweepable device rejected",
                  DCSweep::makeTarget(c, "D1", err) == nullptr
                  && err.find("not sweepable") != std::string::npos);
        ok = ok && (g_failures == 0);
    }

    // ---- 4) divider sweep vs analytic ----
    if (master) std::printf("\n4) Divider sweep vs analytic\n");
    {
        Circuit c("div");
        Node* out = nullptr;
        buildDivider(c, 0.0, out);
        const int idx = nodeIndex(c, out);

        DCSweep sweep(c);
        sweep.setVerbose(false);
        sweep.setSource("V1");
        sweep.setVariation(DCSweepVariation::LIN);
        sweep.setRange(0.0, 10.0, 2.5);  // 0, 2.5, 5, 7.5, 10
        DCSweepResult r = sweep.run(use_mpi);
        checkTrue("sweep converges", r.success);
        ok = ok && r.success;
        if (r.success) {
            checkTrue("5 points", r.points.size() == 5);
            checkTrue("5 points solved", r.stats.points_solved == 5);
            checkTrue("no failures", r.stats.failed_points == 0);
            checkTrue("no retries needed", r.stats.retry_solves == 0);
            bool exact = (r.points.size() == 5);
            for (size_t i = 0; exact && i < r.points.size(); ++i) {
                const double vin = r.points[i].sweep_value;
                const double vout = r.points[i].node_voltages[idx];
                if (std::abs(vout - vin / 2.0) > 1e-6) exact = false;
            }
            checkTrue("Vout == Vin/2 at every point", exact);
            checkTrue("branch currents reported",
                      r.points.size() == 5 && !r.points[0].branch_currents.empty());
            // Source restored to the netlist value after the sweep.
            Device* raw = c.getDevice("V1");
            auto* vsrc = dynamic_cast<VoltageSource*>(raw);
            checkTrue("source found after sweep", vsrc != nullptr);
            if (vsrc) checkClose("source restored", vsrc->waveform().value(0.0), 0.0, 0.0);
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 5) resistor sweep (Xyce `.DC R1 ...`) ----
    if (master) std::printf("\n5) Resistor sweep vs analytic\n");
    {
        Circuit c("rsweep");
        Node* out = nullptr;
        buildDivider(c, 10.0, out, 1e3, 1e3);
        const int idx = nodeIndex(c, out);

        DCSweepSpec spec;
        spec.source = "R1";
        spec.variation = DCSweepVariation::LIN;
        spec.start = 1e3;
        spec.stop  = 3e3;
        spec.step  = 1e3;

        DCSweep sweep(c);
        sweep.setVerbose(false);
        sweep.setSweep(spec);
        DCSweepResult r = sweep.run(use_mpi);
        checkTrue("resistor sweep converges", r.success);
        if (r.success) {
            checkTrue("3 points", r.points.size() == 3);
            bool exact = (r.points.size() == 3);
            for (size_t i = 0; exact && i < r.points.size(); ++i) {
                const double rr = r.points[i].sweep_value;          // R1
                const double expected = 10.0 * 1e3 / (rr + 1e3);   // Vout
                const double vout = r.points[i].node_voltages[idx];
                if (std::abs(vout - expected) > 1e-6) exact = false;
            }
            checkTrue("Vout == 10*R2/(R1+R2) at every point", exact);
            Device* raw = c.getDevice("R1");
            auto* r1 = dynamic_cast<Resistor*>(raw);
            checkTrue("resistance restored after sweep",
                      r1 && std::abs(r1->getEffectiveResistance() - 1e3) < 1e-6);
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 6) nested sweep ----
    if (master) std::printf("\n6) Nested sweep\n");
    {
        // V1 a 0 DC 0 (outer), V2 b 0 DC 0 (inner), R1 a b: ideal sources
        // pin their nodes, so v(a) == V1 and v(b) == V2 analytically.
        Circuit c("nested");
        Node* gnd = c.getGroundNode();
        Node* a = c.addNode(1, "a");
        Node* b = c.addNode(2, "b");
        c.addDevice(make_rcp<VoltageSource>("V1", a, gnd, SourceWaveform::dc(0.0)));
        c.addDevice(make_rcp<VoltageSource>("V2", b, gnd, SourceWaveform::dc(0.0)));
        c.addDevice(make_rcp<Resistor>("R1", a, b, 1e3));
        const int ia = nodeIndex(c, a);
        const int ib = nodeIndex(c, b);

        DCSweepSpec outer;
        outer.source = "V1";
        outer.start = 0.0; outer.stop = 2.0; outer.step = 1.0;
        DCSweepSpec inner;
        inner.source = "V2";
        inner.start = 0.0; inner.stop = 2.0; inner.step = 1.0;

        DCSweep sweep(c);
        sweep.setVerbose(false);
        sweep.setSweep(outer);
        sweep.setInnerSweep(inner);
        DCSweepResult r = sweep.run(use_mpi);
        checkTrue("nested sweep converges", r.success);
        if (r.success) {
            checkTrue("nested flag", r.nested && r.sweep_source2 == "V2");
            checkTrue("3x3 = 9 points", r.points.size() == 9);
            // Expected traversal: outer 0 -> inner 0,1,2; outer 1 -> 0,1,2; ...
            bool exact = (r.points.size() == 9);
            for (size_t i = 0; exact && i < r.points.size(); ++i) {
                const DCSweepPoint& p = r.points[i];
                const double eo = static_cast<double>(i / 3);
                const double ei = static_cast<double>(i % 3);
                if (std::abs(p.outer_value - eo) > 1e-12) exact = false;
                if (std::abs(p.sweep_value - ei) > 1e-12) exact = false;
                if (std::abs(p.node_voltages[ia] - eo) > 1e-6) exact = false;
                if (std::abs(p.node_voltages[ib] - ei) > 1e-6) exact = false;
            }
            checkTrue("v(a)=V1 and v(b)=V2 at every point", exact);
            // Inner loop restarts at its own start value on each outer point.
            bool restart = (r.points.size() == 9
                            && std::abs(r.points[3].sweep_value - 0.0) < 1e-12
                            && std::abs(r.points[3].outer_value - 1.0) < 1e-12
                            && std::abs(r.points[6].sweep_value - 0.0) < 1e-12
                            && std::abs(r.points[6].outer_value - 2.0) < 1e-12);
            checkTrue("inner loop restarts per outer point", restart);
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 7) predictor equivalence ----
    if (master) std::printf("\n7) Linear predictor equivalence\n");
    {
        auto runWith = [&](bool predictor, std::vector<double>& out) {
            Circuit c("pred");
            Node* outNode = nullptr;
            buildDivider(c, 0.0, outNode, 1e3, 3e3);
            const int idx = nodeIndex(c, outNode);
            DCSweep sweep(c);
            sweep.setVerbose(false);
            sweep.setSource("V1");
            sweep.setPredictor(predictor);
            sweep.setRange(0.0, 10.0, 0.5);
            DCSweepResult r = sweep.run(use_mpi);
            out.clear();
            for (const auto& p : r.points) out.push_back(p.node_voltages[idx]);
            return r;
        };
        std::vector<double> with_pred, without_pred;
        const DCSweepResult rp = runWith(true, with_pred);
        const DCSweepResult rn = runWith(false, without_pred);
        checkTrue("both sweeps converge", rp.success && rn.success);
        checkTrue("predictor used", rp.stats.predictor_hits > 0);
        checkTrue("predictor off = no hits", rn.stats.predictor_hits == 0);
        bool same = (with_pred.size() == without_pred.size()
                     && !with_pred.empty());
        for (size_t i = 0; same && i < with_pred.size(); ++i) {
            if (std::abs(with_pred[i] - without_pred[i]) > 1e-9) same = false;
        }
        checkTrue("predictor does not change the solution", same);
        ok = ok && (g_failures == 0);
    }

    // ---- 8) failure ladder (source stepping + NOFAIL) ----
    if (master) std::printf("\n8) Failure ladder\n");
    {
        // V1 and V2 in parallel with different values: no DC solution for
        // most sweep values (singular Jacobian), so the sweep either aborts
        // or, with NOFAIL, keeps going and flags the rows.
        auto buildConflict = [](Circuit& c, Node*& a) {
            Node* gnd = c.getGroundNode();
            a = c.addNode(1, "a");
            c.addDevice(make_rcp<VoltageSource>("V1", a, gnd, SourceWaveform::dc(0.0)));
            c.addDevice(make_rcp<VoltageSource>("V2", a, gnd, SourceWaveform::dc(3.0)));
            c.addDevice(make_rcp<Resistor>("R1", a, gnd, 1e3));
        };
        {
            Circuit c("conflict-abort");
            Node* a = nullptr;
            buildConflict(c, a);
            DCSweep sweep(c);
            sweep.setVerbose(false);
            sweep.setSource("V1");
            sweep.setRange(0.0, 4.0, 1.0);
            DCSweepResult r = sweep.run(use_mpi);
            checkTrue("non-convergence aborts by default", !r.success);
            checkTrue("abort names the failing point",
                      r.message.find("DC sweep failed at V1=") != std::string::npos);
            // The source is restored even on the error path.
            auto* vsrc = dynamic_cast<VoltageSource*>(c.getDevice("V1"));
            checkTrue("source restored after abort",
                      vsrc && std::abs(vsrc->waveform().value(0.0)) < 1e-12);
        }
        {
            Circuit c("conflict-nofail");
            Node* a = nullptr;
            buildConflict(c, a);
            const int ia = nodeIndex(c, a);
            DCSweep sweep(c);
            sweep.setVerbose(false);
            sweep.setSource("V1");
            sweep.setRange(0.0, 4.0, 1.0);
            sweep.setContinueOnFail(true);
            DCSweepResult r = sweep.run(use_mpi);
            checkTrue("NOFAIL keeps the sweep alive", r.success);
            checkTrue("failures counted", r.stats.failed_points > 0);
            checkTrue("source-step refinement exercised", r.stats.retry_solves > 0);
            checkTrue("all points attempted", r.stats.points_solved == 5);
            bool flagged = false;
            for (const auto& p : r.points) if (!p.converged) flagged = true;
            checkTrue("failed rows flagged", flagged);
            checkTrue("NOFAIL reported in message",
                      r.message.find("without solution") != std::string::npos);
            // A row without solution holds the last converged values.
            double last_good = -1.0;
            bool holds = true;
            for (const auto& p : r.points) {
                if (!p.converged) {
                    if (last_good >= 0.0
                        && std::abs(p.node_voltages[ia] - last_good) > 1e-6) {
                        holds = false;
                    }
                } else {
                    last_good = p.node_voltages[ia];
                }
            }
            checkTrue("failed row holds last converged solution", holds);
        }
        ok = ok && (g_failures == 0);
    }

    // ---- 9) streaming CSV + bounded RAM ----
    if (master) std::printf("\n9) Streaming CSV\n");
    {
        Circuit c("csv");
        Node* out = nullptr;
        buildDivider(c, 0.0, out, 1e3, 2e3);   // Vout = Vin/3 -> needs digits
        const std::vector<std::string> labels{"Vin", "Vout"};
        std::ostringstream csv;
        csv.precision(6);   // DCSweep must raise it and restore it

        DCSweep sweep(c);
        sweep.setVerbose(false);
        sweep.setSource("V1");
        sweep.setRange(0.0, 8.0, 0.5);            // 17 points
        sweep.setCsvStream(&csv, labels);
        sweep.setMaxRetainedPoints(5);            // RAM cap
        DCSweepResult r = sweep.run(use_mpi);
        checkTrue("csv sweep converges", r.success);

        std::istringstream in(csv.str());
        std::string header, line;
        std::getline(in, header);
        checkTrue("csv header", header == "V1,Vin,Vout");
        size_t rows = 0;
        size_t fields = 0;
        while (std::getline(in, line)) { ++rows; fields = std::count(line.begin(), line.end(), ',') + 1; }
        checkTrue("csv has every point", rows == 17);
        checkTrue("csv row width", fields == 3);
        checkTrue("csv full precision", csv.str().find("0.333333333") != std::string::npos);
        checkTrue("in-memory sample decimated", r.points.size() <= 6 && !r.points.empty());
        checkTrue("stream precision restored", csv.precision() == 6);
        checkClose("sample starts at 0", r.points.front().sweep_value, 0.0, 0.0);
        checkClose("sample ends at 8", r.points.back().sweep_value, 8.0, 0.0);
        ok = ok && (g_failures == 0);
    }

    // ---- 10) temperature: `.TEMP`, `.DC TEMP ...`, card precedence ----
    if (master) std::printf("\n10) Circuit temperature\n");
    {
        // (a) parser: `.TEMP` / `.OPTIONS TEMP=` / `SWEEPPAR`
        {
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\n.MODEL DMOD D (IS=1e-14)\nV1 in 0 DC 5\nD1 in 0 DMOD\n"
                ".TEMP 85\n.OP\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue(".TEMP parses", parsed && p.temp_given);
            checkClose(".TEMP value (Celsius)", p.temp_value, 85.0, 0.0);
        }
        {
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 in 0 DC 5\nR1 in 0 1k\n"
                ".OPTIONS TEMP=-40C SWEEPPAR\n.OP\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue(".OPTIONS TEMP/SWEEPPAR parse",
                      parsed && p.temp_given && p.opt_sweep_parallel);
            checkClose(".OPTIONS TEMP value", p.temp_value, -40.0, 0.0);
        }
        {
            Circuit c("p");
            NetlistParser parser;
            const bool parsed = parser.parseString(
                "* t\nV1 in 0 DC 5\nR1 in 0 1k\n.TEMP 300.15K\n.OP\n.END\n", c);
            const auto& p = parser.getAnalysisParameters();
            checkTrue("Kelvin marker accepted",
                      parsed && std::abs(p.temp_value - 27.0) < 1e-9);
        }

        // (b) R(T) follows the card's TC1, with TNOM as the reference:
        //     Vout = 10 * R2 / (R1(T) + R2), R1(T) = 1k*(1 + TC1*(T-TNOM))
        auto dividerWithTC = [](Circuit& c, double /*tc1*/, Node*& out,
                                RCP<const ResistorModel> card = RCP<const ResistorModel>()) {
            Node* gnd = c.getGroundNode();
            Node* in = c.addNode(1, "in");
            out = c.addNode(2, "out");
            c.addDevice(make_rcp<VoltageSource>("V1", in, gnd,
                                                SourceWaveform::dc(10.0)));
            c.addDevice(make_rcp<Resistor>("R1", in, out, 1e3, card));
            c.addDevice(make_rcp<Resistor>("R2", out, gnd, 1e3));
        };
        auto expectedVout = [](double tc1, double tC) {
            const double r = 1e3 * (1.0 + tc1 * (tC - 27.0));
            return 10.0 * 1e3 / (r + 1e3);
        };

        {   // (b1) static .TEMP through the circuit
            const double tc1 = 0.0039;
            Circuit c("temp");
            Node* out = nullptr;
            RCP<ResistorModel> card = make_rcp<ResistorModel>("RCU");
            card->TC1 = tc1;
            dividerWithTC(c, tc1, out, card);
            const int idx = nodeIndex(c, out);

            c.setTemperature(273.15 + 85.0);
            checkClose("circuit temperature stored", c.temperature(),
                       273.15 + 85.0, 0.0);
            DCAnalysis dc(c);
            dc.setVerbose(false);
            DCResult r = dc.run(use_mpi, false);
            checkTrue(".TEMP DC converges", r.success);
            if (r.success) {
                checkClose("R(T) divider at 85 C", r.node_voltages[idx],
                           expectedVout(tc1, 85.0), 1e-6);
            }
            ok = ok && (g_failures == 0);
        }

        {   // (b2) a card that pins TEMP= ignores the circuit temperature
            //      (SPICE precedence: instance TEMP > card TEMP > circuit).
            //      Two identical diode circuits: one card follows the
            //      circuit temperature, the other stays at its own 27 C.
            auto diodeCircuit = [](ModelRegistry& reg, const std::string& name,
                                   bool pin_temp, Circuit& c) {
                ModelCard card;
                card.name = name;
                card.type = "D";
                card.params["IS"] = 1e-14;
                card.params["N"]  = 1.0;
                if (pin_temp) card.params["TEMP"] = 27.0;   // Celsius on card
                reg.createModel(card);
                Node* gnd = c.getGroundNode();
                Node* a = c.addNode(1, "a");
                Node* b = c.addNode(2, "b");
                c.addDevice(make_rcp<VoltageSource>("V1", a, gnd,
                                                    SourceWaveform::dc(0.6)));
                c.addDevice(make_rcp<Resistor>("R1", a, b, 1e3));
                c.addDevice(make_rcp<Diode>("D1", b, gnd,
                                            reg.findAs<DiodeModel>(name)));
            };
            auto diodeCurrent = [&](Circuit& c, bool hot) {
                if (hot) c.setTemperature(273.15 + 85.0);
                DCAnalysis dc(c);
                dc.setVerbose(false);
                DCResult r = dc.run(use_mpi, false);
                if (!r.success) return -1.0;
                auto* r1 = dynamic_cast<Resistor*>(c.getDevice("R1"));
                return r1 ? r1->getCurrent() : -1.0;
            };

            ModelRegistry regFollow, regPin;
            Circuit follow("t-follow"), pin("t-pin");
            diodeCircuit(regFollow, "DFOLLOW", false, follow);
            diodeCircuit(regPin, "DPINNED", true, pin);

            const double i_ref_follow = diodeCurrent(follow, false);
            const double i_ref_pin    = diodeCurrent(pin,    false);
            const double i_hot_follow  = diodeCurrent(follow, true);
            const double i_hot_pin     = diodeCurrent(pin,    true);
            checkTrue("diode circuits converge",
                      i_ref_follow > 0.0 && i_ref_pin > 0.0);
            checkTrue("card without TEMP follows the circuit temperature",
                      i_hot_follow > i_ref_follow * 1.05);
            checkClose("card with TEMP= ignores the circuit temperature",
                       i_hot_pin, i_ref_pin, 1e-12);
            ok = ok && (g_failures == 0);
        }

        {   // (b3) TEMP as the .DC sweep variable (Xyce), Celsius values
            const double tc1 = 0.0039;
            Circuit c("tempsweep");
            Node* out = nullptr;
            RCP<ResistorModel> card = make_rcp<ResistorModel>("RCU");
            card->TC1 = tc1;
            dividerWithTC(c, tc1, out, card);
            const int idx = nodeIndex(c, out);

            DCSweepSpec spec;
            spec.source = "TEMP";     // keyword, not a device
            spec.start = -40.0;
            spec.stop  = 85.0;
            spec.step  = 25.0;        // -40, -15, 10, 35, 60, 85 = 6 points

            DCSweep sweep(c);
            sweep.setVerbose(false);
            sweep.setSweep(spec);
            DCSweepResult r = sweep.run(use_mpi);
            checkTrue("TEMP sweep converges", r.success);
            if (r.success) {
                checkTrue("6 temperature points", r.points.size() == 6);
                bool exact = (r.points.size() == 6);
                for (size_t i = 0; exact && i < r.points.size(); ++i) {
                    const DCSweepPoint& p = r.points[i];
                    if (std::abs(p.node_voltages[idx]
                                 - expectedVout(tc1, p.sweep_value)) > 1e-6) {
                        exact = false;
                    }
                }
                checkTrue("Vout == 10*R2/(R1(T)+R2) at every T", exact);
                checkClose("temperature restored after sweep",
                           c.temperature(), 300.15, 0.0);
            }
            ok = ok && (g_failures == 0);
        }

        {   // (b4) an unsweepable name that is not TEMP is still rejected
            Circuit c("bad");
            Node* out = nullptr;
            dividerWithTC(c, 0.0, out);
            std::string err;
            checkTrue("unknown variable still rejected",
                      DCSweep::makeTarget(c, "NOPE", err) == nullptr
                      && err.find("not TEMP") != std::string::npos);
            ok = ok && (g_failures == 0);
        }
    }

    // ---- 11) sweep points distributed across ranks ----
    if (master) std::printf("\n11) Point distribution across ranks\n");
    {
        // Under MPI with .OPTIONS SWEEPPAR each rank solves its own block of
        // the grid, so the answers come from different seeds: the test
        // therefore checks ANALYTIC agreement (the divider is exact) and
        // ordering, not bit-identity. Serial runs must stay bit-reproducible,
        // which is what the other tests already cover.
        Circuit c("dist");
        Node* out = nullptr;
        buildDivider(c, 0.0, out);
        const int idx = nodeIndex(c, out);

        DCSweepSpec spec;
        spec.source = "V1";
        spec.start = 0.0;
        spec.stop  = 10.0;
        spec.step  = 0.5;          // 21 points, >= 2 per rank below

        DCSweep sweep(c);
        sweep.setVerbose(false);
        sweep.setSweep(spec);
        sweep.setPointDistribution(true);
        DCSweepResult r = sweep.run(use_mpi);

        checkTrue("distributed sweep converges", r.success);
        if (r.success) {
            const int ranks = mpi.getSize();
            const int expected_ranks =
                (ranks > 1 && 21 >= 2 * ranks) ? ranks : 1;
            // Statistics are broadcast: every rank reports the same outcome.
            checkTrue("ranks reported on every rank",
                      r.stats.ranks_used == expected_ranks);
            checkTrue("every point solved", r.stats.points_solved == 21);
            // Points and the CSV belong to the master in this mode.
            if (mpi.isMaster()) {
                checkTrue("every point present", r.points.size() == 21);
                bool ordered = true, exact = (r.points.size() == 21);
                for (size_t i = 0; exact && i < r.points.size(); ++i) {
                    const DCSweepPoint& p = r.points[i];
                    if (std::abs(p.sweep_value - 0.5 * i) > 1e-9) ordered = false;
                    if (std::abs(p.node_voltages[idx] - 0.25 * i) > 1e-6) {
                        exact = false;
                    }
                }
                checkTrue("points come back in sweep order", ordered);
                checkTrue("Vout == Vin/2 at every point (all ranks)", exact);
            }
        }
        ok = ok && (g_failures == 0);
    }

    if (master) {
        std::printf(ok && g_failures == 0 ? "DC sweep PASSED\n"
                                          : "DC sweep FAILED (%d)\n", g_failures);
    }
    mpi.finalize();
    return (ok && g_failures == 0) ? 0 : 1;
}
