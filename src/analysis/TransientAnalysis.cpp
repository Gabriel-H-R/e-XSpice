/**
 * @author Gabriel Henrique Silva
 * @file TransientAnalysis.cpp
 * @brief Transient analysis on top of the DAE device layer.
 *
 * Algorithm (Xyce-style):
 *   1. setup(): fixed sparsity pattern, state history, MNA indices;
 *   2. operating point at t = tstart:
 *        - devices with IC= present  → IC operating point (capacitors become
 *          voltage sources, inductors current sources), then the charge
 *          history is frozen from that point;
 *        - otherwise a plain DC OP (with sources at their tstart value);
 *   3. time loop:
 *        for each step dt:
 *          solve the nonlinear DAE at t+dt with
 *              J  = dFdx + s·dQdx,   s = 2/dt (trap) or 1/dt (BE)
 *              RHS = B − F − s·Q + s·Q_n (+ dQ/dt_n for trap)
 *          on failure: halve dt (local step rejection) and retry;
 *          on success: accept states, compute dQ/dt, extract, record;
 *   4. breakpoints (source discontinuities) clamp dt so steps never jump
 *      over a waveform edge — no interpolation artifacts.
 */

#include "TransientAnalysis.h"
#include "DaeSystem.h"
#include "../devices/Capacitor.h"
#include "../devices/Inductor.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>

namespace e_XSpice {

TransientAnalysis::TransientAnalysis(Circuit& circuit)    : circuit_(circuit)
    , method_(IntegrationMethod::TRAPEZOIDAL)
    , solver_type_(SolverType::BICGSTAB)
    , max_iterations_(100)
    , tolerance_(1e-9)
    , verbose_(false)
    , save_interval_(1)
{
    system_ = std::make_unique<DaeSystem>(circuit);
}

TransientAnalysis::~TransientAnalysis() = default;

double TransientAnalysis::estimateLteRatio(double dt,
                                           const std::vector<double>& x_new,
                                           const std::vector<double>& x_old,
                                           const std::vector<double>& x_older,
                                           double reltol, double vntol,
                                           double abstol) const {
    const int n = static_cast<int>(x_new.size());
    double worst = 0.0;
    for (int i = 0; i < n; ++i) {
        const double lte = std::abs(dt / 6.0) *
            std::abs(x_new[i] - 2.0 * x_old[i] + x_older[i]);
        const double tol = reltol * std::abs(x_new[i]) +
                           (system_->isCurrentUnknown(i) ? abstol : vntol);
        if (tol > 0.0) worst = std::max(worst, lte / tol);
    }
    return worst;
}

void TransientAnalysis::createNodeMap(std::vector<int>& node_map) {
    node_map = system_->nodeMap();
}

TransientDataPoint TransientAnalysis::collectDataPoint(double time,
                                                       const std::vector<int>& node_map) {
    (void)node_map;
    TransientDataPoint point;
    point.time = time;
    point.node_voltages.reserve(circuit_.getNodes().size());
    for (Node* n : circuit_.getNodes()) {
        point.node_voltages.push_back(n->getVoltage());
    }
    for (Device* d : circuit_.getDevices()) {
        if (d->numExtraUnknowns() > 0 ||   // V sources, inductors
            d->getType() == DeviceType::CAPACITOR ||
            d->getType() == DeviceType::DIODE) {
            point.branch_currents.push_back(d->getCurrent());
        }
    }
    return point;
}


TransientResult TransientAnalysis::run(double tstart, double tstop, double tstep,
                                       bool use_mpi) {
    (void)use_mpi;
    TransientResult result;
    result.success = false;

    if (tstep <= 0.0) {
        result.message = "Time step must be positive";
        return result;
    }
    if (tstop <= tstart) {
        result.message = "Stop time must be greater than start time";
        return result;
    }

    std::string error_msg;
    if (!circuit_.validate(error_msg)) {
        result.message = "Circuit validation failed: " + error_msg;
        return result;
    }

    if (verbose_) {
        std::cout << "\n===== Transient Analysis (DAE formulation) =====\n";
        std::cout << "Time: " << tstart << "s to " << tstop << "s, step " << tstep << "s\n";
        std::cout << "Method: " << (method_ == IntegrationMethod::TRAPEZOIDAL
                                        ? "Trapezoidal" : "Backward Euler") << "\n";
    }

    // ================= SETUP =================
    system_->setup();
    system_->setMaxIterations(max_iterations_);
    system_->setRelTol(1e-3);
    system_->setAbsTol(std::max(1e-15, tolerance_ * 1e-6));
    system_->setGminSchedule(1e-3, 1e-12);
    system_->setSolverType(solver_type_);
    system_->setIntegrator(method_ == IntegrationMethod::TRAPEZOIDAL);
    system_->setVerbose(verbose_);

    std::vector<int> node_map;
    createNodeMap(node_map);

    // ================= BREAKPOINTS =================
    std::vector<double> breakpoints;
    system_->requestBreakpoints(breakpoints, tstop);
    breakpoints.push_back(tstop);
    std::sort(breakpoints.begin(), breakpoints.end());

    // ================= OPERATING POINT =================
    // Devices with an initial condition force the IC mode of the OP
    // (Xyce dcopFlag + ICGiven convention).
    bool any_ic = false;
    for (const Device* d : circuit_.getDevices()) {
        const DeviceType ty = d->getType();
        if (ty == DeviceType::CAPACITOR || ty == DeviceType::INDUCTOR) {
            if (const auto* c = dynamic_cast<const class Capacitor*>(d)) {
                if (c->hasIC()) any_ic = true;
            }
            if (const auto* i = dynamic_cast<const class Inductor*>(d)) {
                if (i->hasIC()) any_ic = true;
            }
        }
    }

    system_->beginSolve(/*initial_solve=*/true);
    system_->setICOperatingPoint(any_ic);
    {
        std::vector<double> guess(system_->size(), 0.0);
        system_->seedInitialGuess(guess);
        NLSolveResult op = system_->solveNonlinear(tstart, 0.0, false, guess);
        if (!op.converged) {
            result.message = "Operating point failed: " + op.message;
            return result;
        }
        system_->extractSolution();
        // Freeze the charge/flux history from the operating point.
        system_->commitStates();
        system_->initializeHistory();
    }
    system_->setICOperatingPoint(false);

    // ================= TIME LOOP =================
    auto t0 = std::chrono::high_resolution_clock::now();

    result.data.reserve(static_cast<size_t>((tstop - tstart) / tstep) + 2);
    result.data.push_back(collectDataPoint(tstart, node_map));

    double t = tstart;
    int step = 0;
    const int num_steps_estimate =
        static_cast<int>(std::ceil((tstop - tstart) / tstep));

    // Solution history for the LTE estimate (trapezoidal, order 2):
    //   LTE ≈ (h/6)·|x_{n+1} − 2·x_n + x_{n−1}|
    std::vector<double> x_prev_prev = system_->solution();   // t_{n-1}
    std::vector<double> x_prev      = system_->solution();   // t_n

    // Tolerances (SPICE-style: separate voltage/current absolutes)
    const double lte_reltol = 1e-3;
    const double lte_vntol  = std::max(1e-9, tolerance_ * 1e-3);
    const double lte_abstol = std::max(1e-15, tolerance_ * 1e-6);
    const double max_growth = 4.0;

    double dt_want = tstep;

    while (t < tstop - 1e-18 * std::max(1.0, std::abs(tstop))) {
        double dt = std::min(dt_want, tstop - t);

        // Never step over a breakpoint.
        for (double bp : breakpoints) {
            if (bp > t + 1e-15 * std::abs(bp + 1.0) && bp < t + dt) {
                dt = bp - t;
                break;
            }
        }
        // Device-requested maximum step (0.1·PER for PULSE, 0.1/FREQ…)
        const double dev_hint = system_->maxTimeStepHint(t);
        if (dev_hint > 0.0) dt = std::min(dt, dev_hint);

        if (dt <= 0.0) break;

        // ---- try the step (reduce dt on Newton failure or LTE reject) ----
        double ratio = 1.0;
        bool accepted_step = false;
        for (int attempt = 0; attempt < 20; ++attempt) {
            const NLSolveResult nl =
                system_->solveNonlinear(t + dt, integrationCoeff(dt), true, x_prev);

            if (!nl.converged) {
                // Too big (or too hard): shrink and retry from the same
                // accepted point.
                system_->restoreAcceptedState();
                dt *= 0.25;
                if (dt < 1e-15 * std::max(1.0, std::abs(tstop))) break;
                continue;
            }

            // LTE check (skipped while the history is too short)
            if (lte_ready_) {
                ratio = estimateLteRatio(dt, system_->solution(),
                                         x_prev, x_prev_prev,
                                         lte_reltol, lte_vntol, lte_abstol);
            } else {
                ratio = 0.0;
            }

            if (ratio > 1.0) {
                system_->restoreAcceptedState();
                dt *= std::max(0.1, 0.9 * std::cbrt(1.0 / ratio));
                continue;
            }
            accepted_step = true;
            break;
        }

        if (!accepted_step) {
            result.message = "Time step failed at t=" + std::to_string(t) + "s";
            return result;
        }

        // ---- accept the step ----
        t += dt;
        ++step;
        if (step >= 2) lte_ready_ = true;
        if (lte_ready_) {
            x_prev_prev = x_prev;
        }
        x_prev = system_->solution();

        // Advance the charge/flux history with the integrator formula.
        system_->acceptTimeStep(integrationCoeff(dt),
                                method_ == IntegrationMethod::TRAPEZOIDAL);
        system_->extractSolution();

        // Next step: follow the LTE ratio (bounded growth).
        const double growth = (ratio > 0.0)
            ? std::min(max_growth, std::max(0.5, 0.9 * std::cbrt(1.0 / ratio)))
            : max_growth;
        dt_want = std::min(tstep, dt * growth);

        if (step % save_interval_ == 0) {
            result.data.push_back(collectDataPoint(t, node_map));
        }
        if (progress_callback_) {
            progress_callback_(t, step, num_steps_estimate);
        }
        if (verbose_ && step % 200 == 0) {
            std::cout << "  step " << step << " (t=" << t << "s, dt=" << dt << ")\n";
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    result.total_time = std::chrono::duration<double>(t1 - t0).count();
    result.success = true;
    result.message = "Transient analysis completed successfully";

    if (verbose_) {
        std::cout << "Completed " << step << " steps in " << result.total_time << " s\n";
    }
    return result;
}

} // namespace e_XSpice
