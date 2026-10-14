/**
 * @author Gabriel Henrique Silva
 * @file DCAnalysis.cpp
 * @brief DC operating point analysis on top of the DAE device layer.
 */

#include "DCAnalysis.h"
#include "DaeSystem.h"
#include <iostream>
#include <chrono>
#include <algorithm>

namespace e_XSpice {

DCAnalysis::DCAnalysis(Circuit& circuit)
    : circuit_(circuit)
    , solver_type_(SolverType::BICGSTAB)
    , max_iterations_(200)
    , tolerance_(1e-9)
    , verbose_(false)
{
}

DCResult DCAnalysis::run(bool use_mpi, bool use_initial_guess) {
    DCResult result;
    (void)use_mpi;         // DC is a single dense-reduction ep; OpenMP is
    (void)use_initial_guess; // used inside the solver. MPI: future work.

    std::string error_msg;
    if (!circuit_.validate(error_msg)) {
        result.success = false;
        result.message = "Circuit validation failed: " + error_msg;
        return result;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    DaeSystem system(circuit_);
    system.setup();
    system.setMaxIterations(max_iterations_);
    system.setRelTol(1e-3);
    system.setAbsTol(std::max(1e-12, tolerance_ * 1e-3));
    system.setGminSchedule(1e-3, 1e-12);
    system.setSolverType(solver_type_);
    system.setVerbose(verbose_);
    system.beginSolve(/*initial_solve=*/true);

    // Sources evaluated at t = 0 (DC value of every waveform).
    const NLSolveResult nl = system.solveNonlinear(/*time=*/0.0,
                                                   /*q_coeff=*/0.0,
                                                   /*use_history=*/false);

    result.assembly_time = nl.elapsed_seconds;
    result.solve_time = 0.0;

    if (!nl.converged) {
        result.success = false;
        result.message = "Newton did not converge: " + nl.message;
        return result;
    }

    result.solver_result.converged = true;
    result.solver_result.iterations = nl.total_iterations;
    result.solver_result.residual_norm = nl.residual_norm;
    result.solver_result.elapsed_time = nl.elapsed_seconds;
    result.solver_result.message = nl.message;

    system.extractSolution();

    // ---- collect outputs ----
    auto t1 = std::chrono::high_resolution_clock::now();
    result.assembly_time = 0.0;
    result.solve_time = std::chrono::duration<double>(t1 - t0).count();

    result.node_voltages.clear();
    result.branch_currents.clear();
    for (Node* n : circuit_.getNodes()) {
        result.node_voltages.push_back(n->getVoltage());
    }
    for (Device* d : circuit_.getDevices()) {
        result.branch_currents.push_back(d->getCurrent());
    }

    result.success = true;
    result.message = "DC operating point converged";

    if (verbose_) {
        std::cout << "[DC] converged in " << nl.total_iterations
                  << " iterations (" << nl.gmin_steps << " gmin steps), "
                  << result.solve_time * 1e3 << " ms\n";
    }

    return result;
}

} // namespace e_XSpice
