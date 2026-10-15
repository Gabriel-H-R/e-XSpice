/**
 * @file test_mpi_ladder.cpp
 * @brief Serial-vs-distributed equivalence on an RC ladder (MPI + OpenMP).
 *
 * Builds a 150-section RC ladder driven by a PULSE source (301 devices,
 * ~152 unknowns — above the distribution thresholds) and runs the SAME
 * transient twice in one process:
 *   1. EXSPICE_DISTRIBUTE=0  → fully serial reference;
 *   2. EXSPICE_DISTRIBUTE=1  → device-partitioned assembly + distributed
 *      Krylov solve (active when the job runs under mpiexec -n 2+).
 *
 * The final node voltages must match to solver tolerance. Under plain
 * ctest (size==1) both runs take the serial path and the test trivially
 * passes; under `mpiexec -n 2` it validates the whole distributed stack:
 * assembly Allreduces, state-history sync, convergence votes, distributed
 * BiCGSTAB, LTE/dt control-flow identity and output refresh.
 */

#include "../src/core/Circuit.h"
#include "../src/core/Node.h"
#include "../src/devices/Resistor.h"
#include "../src/devices/VoltageSource.h"
#include "../src/devices/Capacitor.h"
#include "../src/devices/SourceWaveform.h"
#include "../src/analysis/TransientAnalysis.h"
#include "../src/parallel/MPIManager.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace e_XSpice;

namespace {

void setDistEnv(const char* value) {
#ifdef _WIN32
    _putenv_s("EXSPICE_DISTRIBUTE", value);
#else
    setenv("EXSPICE_DISTRIBUTE", value, 1);
#endif
}

constexpr int kSections = 150;

// Builds the ladder: V1(1) -R- node -R- ... with C shunts to ground.
void buildLadder(Circuit& c) {
    Node* gnd = c.getGroundNode();
    Node* prev = c.addNode(1, "1");
    c.addDevice(make_rcp<VoltageSource>(
        "V1", prev, gnd,
        SourceWaveform::pulse(0.0, 5.0, 0.0, 10e-6, 10e-6, 500e-6, 1e-3)));
    for (int i = 1; i <= kSections; ++i) {
        Node* next = c.addNode(i + 1, std::to_string(i + 1));
        c.addDevice(make_rcp<Resistor>("R" + std::to_string(i), prev, next, 1000.0));
        c.addDevice(make_rcp<Capacitor>("C" + std::to_string(i), next, gnd, 10e-9,
                                  0.0, false));
        prev = next;
    }
}

bool runLadder(std::vector<double>& final_voltages, int& num_points) {
    Circuit circuit("Ladder");
    buildLadder(circuit);

    TransientAnalysis tran(circuit);
    tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
    tran.setTolerance(1e-9);
    tran.setSaveInterval(20);
    tran.setVerbose(false);

    TransientResult r = tran.run(0.0, 1e-3, 5e-6, true);
    if (!r.success || r.data.empty()) {
        std::printf("run failed: %s\n", r.message.c_str());
        return false;
    }
    num_points = static_cast<int>(r.data.size());
    final_voltages = r.data.back().node_voltages;
    return true;
}

} // namespace

int main(int argc, char** argv) {
    MPIManager& mpi = MPIManager::getInstance();
    mpi.initialize(&argc, &argv);

    // "auto" mode: single run with the ambient EXSPICE_DISTRIBUTE setting
    // (used by ctest to exercise the DEFAULT thresholds under mpiexec).
    if (argc > 1 && std::string(argv[1]) == "auto") {
        Circuit circuit("Ladder");
        buildLadder(circuit);
        TransientAnalysis tran(circuit);
        tran.setIntegrationMethod(IntegrationMethod::TRAPEZOIDAL);
        tran.setTolerance(1e-9);
        tran.setSaveInterval(20);
        tran.setVerbose(false);
        TransientResult r = tran.run(0.0, 1e-3, 5e-6, true);
        bool ok = r.success && !r.data.empty();
        if (ok) {
            for (double v : r.data.back().node_voltages) {
                if (!std::isfinite(v) || std::abs(v) > 1e6) { ok = false; break; }
            }
        }
        if (mpi.isMaster()) {
            std::printf(ok ? "MPI ladder auto PASSED\n"
                           : "MPI ladder auto FAILED: %s\n",
                        ok ? "" : r.message.c_str());
        }
        mpi.finalize();
        return ok ? 0 : 1;
    }

    if (mpi.isMaster()) {
        std::printf("MPI ladder equivalence: size=%d\n", mpi.getSize());
    }

    setDistEnv("0");
    std::vector<double> v_serial;
    int n_serial = 0;
    if (!runLadder(v_serial, n_serial)) {
        std::printf("FAIL: serial reference run\n");
        mpi.finalize();
        return 1;
    }

    setDistEnv("1");
    std::vector<double> v_dist;
    int n_dist = 0;
    if (!runLadder(v_dist, n_dist)) {
        std::printf("FAIL: distributed run\n");
        mpi.finalize();
        return 1;
    }
    setDistEnv("0");

    bool ok = true;
    if (v_serial.size() != v_dist.size()) {
        std::printf("FAIL: size mismatch %llu vs %llu\n",
                    (unsigned long long)v_serial.size(),
                    (unsigned long long)v_dist.size());
        ok = false;
    } else {
        double worst = 0.0;
        for (size_t i = 0; i < v_serial.size(); ++i) {
            const double tol = 1e-6 * std::abs(v_serial[i]) + 1e-9;
            const double err = std::abs(v_serial[i] - v_dist[i]);
            worst = std::max(worst, err / tol);
            if (err > tol) {
                std::printf("FAIL: node %llu serial=%.10g dist=%.10g\n",
                            (unsigned long long)i, v_serial[i], v_dist[i]);
                ok = false;
                break;
            }
        }
        if (mpi.isMaster()) {
            std::printf("worst normalized deviation: %.3g (points %d vs %d)\n",
                        worst, n_serial, n_dist);
        }
        if (n_serial != n_dist && mpi.isMaster()) {
            std::printf("note: saved-point counts differ (%d vs %d)\n",
                        n_serial, n_dist);
        }
    }

    if (mpi.isMaster()) {
        std::printf(ok ? "MPI ladder equivalence PASSED\n"
                       : "MPI ladder equivalence FAILED\n");
    }
    mpi.finalize();
    return ok ? 0 : 1;
}
