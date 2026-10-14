/**
 * @author Gabriel Henrique Silva
 * @file TransientAnalysis.h
 * @brief Transient analysis (temporal simulation)
 * 
 * Simulates the circuit behavior over time, solving:
 * G(t) × v(t) = b(t)
 * 
 * For each time step t:
 * 1 - Updates companion models (C, L → equivalent resistors + sources)
 * 2 - Assembles modified MNA system
 * 3 - Solves linear system
 * 4 - Updates device states
 * 
 * Integration methods:
 * - Backward Euler: Implicit, stable, but with truncation error
 * - Trapezoidal: Implicit, more accurate, A-stable
 * 
 * Parallelization:
 * - OpenMP: Parallelizes assembly and solver in each step
 * - MPI: Distributes devices among processes
 */

#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include "../core/Circuit.h"
#include "../solver/SparseMatrix.h"
#include "../solver/LinearSolver.h"
#include "../parallel/MPIManager.h"
#include <vector>
#include <string>
#include <functional>
#include <memory>

namespace e_XSpice {

class DaeSystem;   // defined in DaeSystem.h

/**
 * @brief Transient simulation data point
 */
struct TransientDataPoint {
    double time;                          // Time (s)
    std::vector<double> node_voltages;    // Node voltages
    std::vector<double> branch_currents;  // Currents (V sources, inductors)
};

/**
 * @brief Transient analysis result
 */
struct TransientResult {
    bool success;
    std::string message;
    double total_time;                    // Total simulation time (s)
    std::vector<TransientDataPoint> data; // Temporal history
};
/**
 * @brief Temporal integration method
 */
enum class IntegrationMethod {
    BACKWARD_EULER,  // θ = 1 (fully implicit)
    TRAPEZOIDAL      // θ = 0.5 (Crank-Nicolson)
};

/**
 * @brief Transient analysis
 */
class TransientAnalysis {
public:
    /**
     * @brief Constructor
     * @param circuit Reference to the circuit
     */
    TransientAnalysis(Circuit& circuit);
    ~TransientAnalysis();   // out-of-line (DaeSystem is incomplete here)
    
    /**
     * @brief Runs transient analysis
     * @param tstart Initial time (s)
     * @param tstop Final time (s)
     * @param tstep Time step (s)
     * @param use_mpi Use MPI parallelization
     * @return Analysis result
     */
    TransientResult run(double tstart, double tstop, double tstep,
                       bool use_mpi = true);
    
    /**
     * @brief Defines callback for real-time monitoring
     * Called at each time step
     */
    void setProgressCallback(std::function<void(double time, int step, int total_steps)> callback) {
        progress_callback_ = callback;
    }
    
    /**
     * @brief Defines data saving interval
     * By default, saves all steps. With save_interval > 1, saves memory.
     */
    void setSaveInterval(int interval) { save_interval_ = interval; }
    
    // Configurations
    void setIntegrationMethod(IntegrationMethod method) { method_ = method; }
    void setSolverType(SolverType type) { solver_type_ = type; }
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
    Circuit& circuit_;
    IntegrationMethod method_;
    SolverType solver_type_;
    int max_iterations_;
    double tolerance_;
    bool verbose_;
    int save_interval_;
    
    std::function<void(double, int, int)> progress_callback_;

    std::unique_ptr<DaeSystem> system_;   // DAE/NL engine (owns the system)

    /**
     * @brief Node mapping of the DAE system (node id -> matrix row).
     */
    void createNodeMap(std::vector<int>& node_map);

    /**
     * @brief Collects a data point (node voltages + branch currents) at a
     *        given time.
     */
    TransientDataPoint collectDataPoint(double time,
                                       const std::vector<int>& node_map);

    /**
     * @brief Integration coefficient s applied by the effective Jacobian:
     *        2/dt (trapezoidal) or 1/dt (backward Euler).
     */
    double integrationCoeff(double dt) const {
        return (method_ == IntegrationMethod::TRAPEZOIDAL) ? (2.0 / dt) : (1.0 / dt);
    }

    /**
     * @brief Trapezoidal local truncation error estimate, normalized:
     *   ratio = max_i |LTE_i| / (reltol·|x_i| + vtol_or_atol)
     * A ratio ≤ 1 accepts the step. Uses the last three accepted solutions.
     */
    double estimateLteRatio(double dt,
                            const std::vector<double>& x_new,
                            const std::vector<double>& x_old,
                            const std::vector<double>& x_older,
                            double reltol, double vntol, double abstol) const;

    bool lte_ready_ = false;   // true once three solutions are available
};

} // namespace e_XSpice

#endif // TRANSIENT_ANALYSIS_H
