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

namespace e_XSpice {

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
    
    ~TransientAnalysis() = default;
    
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
    
    /**
     * @brief Executes one time step
     * @param t Current time
     * @param dt Time step
     * @param node_map Node mapping
     * @param use_mpi Use MPI
     * @return true if success
     */
    bool timeStep(double t, double dt, const std::vector<int>& node_map, bool use_mpi);
    
    /**
     * @brief Assembles MNA system for transient analysis
     * Uses companion models for capacitors and inductors
     */
    void assembleTransientMNA(SparseMatrix& G,
                             std::vector<double>& b,
                             const std::vector<int>& node_map,
                             double timestep,
                             bool use_mpi);
    
    /**
     * @brief Creates node mapping and assigns indices
     */
    void createNodeMap(std::vector<int>& node_map);
    
    /**
     * @brief Extracts solution and updates nodes
     */
    void extractSolution(const std::vector<double>& x,
                         const std::vector<int>& node_map);
    
    /**
     * @brief Updates states of all reactive devices
     */
    void updateDeviceStates();
    
    /**
     * @brief Collects data from current step
     */
    TransientDataPoint collectDataPoint(double time, const std::vector<int>& node_map);
};

} // namespace e_XSpice

#endif // TRANSIENT_ANALYSIS_H
