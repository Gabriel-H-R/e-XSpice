/**
 * @author Gabriel Henrique Silva
 * @file LinearSolver.h
 * @brief Linear solvers for sparse systems Ax = b
 * 
 * Implementations:
 * 1. Conjugate Gradient (CG): For symmetric positive-definite matrices
 * 2. BiCGSTAB: For non-symmetric matrices (MNA with voltage sources)
 * 
 * OpenMP parallelization for vector operations and matrix-vector product.
 */

#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "SparseMatrix.h"
#include <vector>
#include <string>

namespace e_XSpice {

/**
 * @brief Solver type
 */
enum class SolverType {
    CONJUGATE_GRADIENT,  // For symmetric matrices
    BICGSTAB,           // For non-symmetric matrices
    DIRECT_LU           // LU decomposition (future)
};

/**
 * @brief Solver result
 */
struct SolverResult {
    bool converged;           // Converged?
    int iterations;          // Number of iterations
    double residual_norm;    // Final residual norm
    double elapsed_time;     // Execution time (s)
    std::string message;     // Status message
};

/**
 * @brief Linear solver class
 */
class LinearSolver {
public:
    /**
     * @brief Constructor
     * @param type Solver type
     */
    LinearSolver(SolverType type = SolverType::BICGSTAB);
    
    ~LinearSolver() = default;
    
    /**
     * @brief Solves system Ax = b
     * @param A Sparse matrix (must be in CSR format)
     * @param b Right-hand side vector
     * @param x Solution vector (output)
     * @param x0 Initial solution (optional)
     * @return Solver result
     */
    SolverResult solve(const SparseMatrix& A,
                      const std::vector<double>& b,
                      std::vector<double>& x,
                      const std::vector<double>* x0 = nullptr);
    
    // Settings
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool verbose) { verbose_ = verbose; }
    
    int getMaxIterations() const { return max_iterations_; }
    double getTolerance() const { return tolerance_; }

private:
    SolverType type_;
    int max_iterations_;
    double tolerance_;
    bool verbose_;
    
    /**
     * @brief Conjugate Gradient
     * Works for symmetric positive-definite matrices
     */
    SolverResult solveConjugateGradient(const SparseMatrix& A,
                                        const std::vector<double>& b,
                                        std::vector<double>& x,
                                        const std::vector<double>* x0);
    
    /**
     * @brief BiCGSTAB (Bi-Conjugate Gradient Stabilized)
     * Works for non-symmetric matrices (general MNA case)
     */
    SolverResult solveBiCGSTAB(const SparseMatrix& A,
                               const std::vector<double>& b,
                               std::vector<double>& x,
                               const std::vector<double>* x0);
    
    // Parallelized vector operations
    double dotProduct(const std::vector<double>& a,
                     const std::vector<double>& b) const;
    
    double vectorNorm(const std::vector<double>& v) const;
    
    void vectorScale(std::vector<double>& v, double alpha) const;
    
    void vectorAdd(std::vector<double>& result,
                  const std::vector<double>& a,
                  double alpha,
                  const std::vector<double>& b) const;
};

} // namespace e_XSpice

#endif // LINEAR_SOLVER_H
