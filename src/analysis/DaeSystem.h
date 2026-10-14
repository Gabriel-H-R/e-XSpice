/**
 * @author Gabriel Henrique Silva
 * @file DaeSystem.h
 * @brief Shared DAE/NL solver core for all analysis engines.
 *
 * Responsibilities
 * ----------------
 *  SETUP (once per analysis):
 *   - assigns the MNA layout: node rows + one extra unknown per
 *     voltage-defined device (voltage sources, inductors);
 *   - lets every device register its charge/flux states in a StateStore;
 *   - builds the FIXED sparsity pattern (union of the dFdx and dQdx
 *     stamps) and hands every device its own NNZ slots;
 *
 *  SOLVE (per Newton-Raphson call):
 *   - assembles F, Q, B, dFdx, dQdx through the devices (no dynamic_cast,
 *     no per-iteration trees/maps — values are written into preallocated
 *     CSR slots, and the Jacobian assembly is parallelizable);
 *   - forms the effective system J = dFdx + s·dQdx,
 *     RHS = B − F − s·(Q − Q_history) (+ dQ/dt history for trapezoidal);
 *   - runs the Newton loop with SPICE-style convergence tests and GMIN
 *     stepping as a homotopy fallback.
 *
 * This is where the "companion model" coefficients come from: the devices
 * never see them (Xyce architecture).
 */

#ifndef DAE_SYSTEM_H
#define DAE_SYSTEM_H

#include "../core/Circuit.h"
#include "../core/DaeSolution.h"
#include "../solver/SparseMatrix.h"
#include "../solver/DaeMatrix.h"
#include "../solver/StateStore.h"
#include "../solver/LinearSolver.h"
#include <vector>
#include <string>

namespace e_XSpice {

/**
 * @brief Outcome of a nonlinear solve.
 */
struct NLSolveResult {
    bool        converged = false;
    int         iterations = 0;      // iterations of the last attempt
    int         total_iterations = 0;// including homotopy steps
    int         gmin_steps = 0;      // homotopy attempts used
    double      residual_norm = 0.0;
    double      elapsed_seconds = 0.0;
    std::string message;
};

class DaeSystem {
public:
    explicit DaeSystem(Circuit& circuit);

    /**
     * @brief Assigns indices, registers states, builds the fixed pattern.
     *        Must be called exactly once before any solve.
     */
    void setup();

    /**
     * @brief Total number of unknowns (node rows + extra unknowns).
     */
    int size() const { return n_; }

    /**
     * @brief Node id -> matrix row (-1 for ground).
     */
    const std::vector<int>& nodeMap() const { return node_map_; }

    /**
     * @brief State (charge/flux) history store owned by the devices.
     */
    StateStore& states() { return store_; }
    const StateStore& states() const { return store_; }

    /**
     * @brief Snapshot / restore of the accepted point (used when a time
     *        step is rejected by LTE control or a failed Newton solve).
     */
    void saveAcceptedState();
    void restoreAcceptedState();

    /**
     * @brief Last used solution vector (unknowns).
     */
    const std::vector<double>& solution() const { return x_; }

    /**
     * @brief Runs the nonlinear solve at time t with charge coefficient s
     *        (s = 0 for DC; s = 2/dt or 1/dt in transient).
     *
     * @param update_history if true (transient), the RHS includes the charge
     *        history terms and dQ/dt history of the previous accepted step.
     * @param x0 initial guess (leave empty to reuse the last solution).
     */
    NLSolveResult solveNonlinear(double time,
                                 double q_coeff,
                                 bool use_history,
                                 const std::vector<double>& x0 = {});

    /**
     * @brief Stores the current Q(x) of every device into the "next" state
     *        slot (called internally after each Newton iteration).
     */
    void commitStates();

    /**
     * @brief Propagates an accepted step: next -> current history, and
     *        computes dQ/dt per state slot with the integrator formula.
     *        Must be called exactly once per accepted time step.
     *
     * @param q_coeff  the s used in the solve (1/dt or 2/dt)
     * @param trapezoidal  selects the trapezoidal dQdt reconstruction
     */
    void acceptTimeStep(double q_coeff, bool trapezoidal);

    /**
     * @brief After the operating-point solve: freeze the state history
     *        (accepted charge = current charge, dQ/dt = 0).
     */
    void initializeHistory() { store_.initializeHistory(); }

    /**
     * @brief Writes the solution into node voltages and device outputs.
     */
    void extractSolution();

    /**
     * @brief Marks the beginning of a nonlinear solve (initial_solve=true
     *        only for the very first solve of the simulation).
     */
    void beginSolve(bool initial_solve);

    /**
     * @brief Seeds the initial guess from the devices (IC= handling).
     */
    void seedInitialGuess(std::vector<double>& x) const;

    /**
     * @brief Enables/disables the initial-condition operating point mode
     *        (capacitors become voltage sources, inductors current sources).
     */
    void setICOperatingPoint(bool on);

    /**
     * @brief Collects the transient breakpoints requested by devices
     *        (source waveform discontinuities up to tstop).
     */
    void requestBreakpoints(std::vector<double>& times, double tstop) const;

    /**
     * @brief Smallest device-requested maximum time step (0 = none).
     */
    double maxTimeStepHint(double t) const;

    /**
     * @brief Sets number of iterations / tolerances.
     */
    void setMaxIterations(int iters) { max_iterations_ = iters; }
    void setRelTol(double tol) { reltol_ = tol; }
    void setAbsTol(double tol) { abstol_ = tol; }
    void setGminSchedule(double gmin_start, double gmin_final) {
        gmin_start_ = gmin_start; gmin_final_ = gmin_final;
    }
    void setSolverType(SolverType t) { solver_type_ = t; }
    void setIntegrator(bool trapezoidal) { trapezoidal_ = trapezoidal; }
    void setVerbose(bool v) { verbose_ = v; }
    void setIterationTrace(bool v) { very_verbose_ = v; }

    /**
     * @brief Total branch-current slots (extra unknowns).
     */
    int numExtraUnknowns() const { return num_extra_; }

    /**
     * @brief Number of node-voltage rows (currents occupy the rest).
     */
    int numNodeRows() const { return num_nodes_; }

    /**
     * @brief True when the unknown at `row` is a branch current (needs the
     *        ABSTOL-style current tolerance in LTE checks).
     */
    bool isCurrentUnknown(int row) const { return row >= num_nodes_; }

    /**
     * @brief Debug access to the last assembled Jacobian (dFdx).
     */
    const DaeMatrix& jacobianF() const { return jac_f_; }
    const std::vector<double>& residualF() const { return F_; }
    const std::vector<double>& rhs() const { return rhs_; }

private:
    Circuit& circuit_;
    std::vector<Device*> devices_;

    int n_ = 0;                 // total unknowns
    int num_nodes_ = 0;         // node rows
    int num_extra_ = 0;         // extra unknowns (branch currents)

    std::vector<int> node_map_; // node id -> row

    DaeMatrix  jac_f_;           // dFdx
    DaeMatrix  jac_q_;           // dQdx
    std::vector<double> F_, Q_, B_;
    std::vector<double> rhs_;
    std::vector<double> x_;
    std::vector<double> x_prev_; // last iteration's solution (convergence)
    std::vector<double> x_saved_; // accepted-point snapshot (LTE/retry)

    StateStore store_;

    std::vector<int> diagonal_slots_;  // diagonal slot of every row
    double node_gmin_ = 0.0;          // homotopy conductance

    // configuration
    int    max_iterations_ = 100;
    double reltol_ = 1e-3;
    double abstol_ = 1e-9;
    double gmin_start_ = 1e-3;
    double gmin_final_ = 1e-12;
    bool   trapezoidal_ = true;
    bool   verbose_ = false;
    bool   very_verbose_ = false;   // per-iteration trace (diagnostics)

    SolverType solver_type_ = SolverType::BICGSTAB;

    // helpers
    void assembleVectors(double time, double q_coeff);
    void assembleJacobians(double time, double q_coeff);
    void buildRHS(double q_coeff, bool use_history);
    bool checkConvergence(const std::vector<double>& x_new) const;
    bool solveLinearSystem(std::vector<double>& dx);
    void applyGminResidual();
    void applyGminDiagonal();
    void resolveDiagonalSlots();
    double currentTime_ = 0.0;
};

} // namespace e_XSpice

#endif // DAE_SYSTEM_H
