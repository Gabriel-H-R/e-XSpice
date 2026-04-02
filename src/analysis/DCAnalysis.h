/**
 * @author Gabriel Henrique Silva
 * @file DCAnalysis.h
 * @brief DC Operating Point Analysis with MPI + OpenMP
 * 
 * This class implements DC analysis using Modified Nodal Analysis (MNA):
 * 1. Assembles conductance matrix G and vector b
 * 2. Solves linear system Gv = b
 * 3. Extracts voltages at nodes and currents in sources
 * 
 * Parallelization:
 * - MPI: Distributes matrix assembly among processes
 * - OpenMP: Parallelizes assembly loops and matrix-vector product
 */

#ifndef DC_ANALYSIS_H
#define DC_ANALYSIS_H

#include "../core/Circuit.h"
#include "../solver/SparseMatrix.h"
#include "../solver/LinearSolver.h"
#include "../parallel/MPIManager.h"
#include <vector>
#include <string>

namespace e_XSpice {

/**
 * @brief DC analysis result
 */
struct DCResult {
    bool success;                      // Analysis success
    std::string message;               // Status message
    double assembly_time;              // Assembly time (s)
    double solve_time;                 // Solve time (s)
    SolverResult solver_result;        // Solver details
    std::vector<double> node_voltages; // Node voltages
    std::vector<double> branch_currents; // Source currents
};

/**
 * @brief DC Operating Point Analysis
 */
class DCAnalysis {
public:
    /**
     * @brief Constructor
     * @param circuit Reference to the circuit
     */
    DCAnalysis(Circuit& circuit);
    
    ~DCAnalysis() = default;
    
    /**
     * @brief Runs DC Operating Point Analysis
     * @param use_mpi Use MPI parallelization
     * @param use_initial_guess Use previous solution as initial guess
     * @return Analysis result
     */
    DCResult run(bool use_mpi = true, bool use_initial_guess = false);
    
    /**
     * @brief Configures linear solver
     */
    void setSolverType(SolverType type) { solver_type_ = type; }
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
    Circuit& circuit_;
    SolverType solver_type_;
    int max_iterations_;
    double tolerance_;
    bool verbose_;
    
    std::vector<double> previous_solution_;  // For initial guess
    
    /**
     * @brief Assembles complete MNA system
     * @param G Conductance matrix (output)
     * @param b Right-hand side vector (output)
     * @param node_map Mapping of node IDs to indices (output)
     * @param use_mpi If true, uses MPI parallelization
     */
    void assembleMNA(SparseMatrix& G, 
                     std::vector<double>& b,
                     std::vector<int>& node_map,
                     bool use_mpi);
    
    /**
     * @brief Serial version of assembly
     */
    void assembleMNA_Serial(SparseMatrix& G,
                            std::vector<double>& b,
                            std::vector<int>& node_map);
    
    /**
     * @brief MPI parallel version of assembly
     */
    void assembleMNA_MPI(SparseMatrix& G,
                         std::vector<double>& b,
                         std::vector<int>& node_map);
    
    /**
     * @brief Extracts solution from vector x and updates nodes
     */
    void extractSolution(const std::vector<double>& x,
                         const std::vector<int>& node_map);
    
    /**
     * @brief Creates mapping of nodes to matrix indices
     * Ground (id=0) does not enter the matrix, so it is mapped to -1
     */
    void createNodeMap(std::vector<int>& node_map);
    
    /**
     * @brief Assigns indices for voltage source currents
     */
    void assignVoltageSourceIndices(std::vector<int>& node_map);
};

} // namespace e_XSpice

#endif // DC_ANALYSIS_H
