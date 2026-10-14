/**
 * @author Gabriel Henrote Silva
 * @file MatrixOps.h
 * @brief Minimal matrix interface consumed by the linear solvers.
 *
 * The Krylov solvers (CG / BiCGSTAB) only need the size and the sparse
 * matrix-vector product. Exposing that single-operation interface lets
 * both the legacy SparseMatrix (map-based) and the new fixed-pattern
 * DaeMatrix be solved by the same code without copies or adapters.
 */

#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <vector>

namespace e_XSpice {

class MatrixOps {
public:
    virtual ~MatrixOps() = default;

    /**
     * @brief Number of rows (systems are always square).
     */
    virtual int rows() const = 0;

    /**
     * @brief y = A·x (parallelized inside the implementations).
     */
    virtual void multiplyVector(const std::vector<double>& x,
                                std::vector<double>& y) const = 0;

    /**
     * @brief Reads a single entry (0.0 when absent). Needed by the direct
     *        dense fallback solver and by tests/debug output.
     */
    virtual double get(int row, int col) const = 0;

    /**
     * @brief True when the matrix is symmetric positive definite
     *        (permits Conjugate Gradient). Default: unknown/false.
     */
    virtual bool isSymmetric() const { return false; }
};

} // namespace e_XSpice

#endif // MATRIX_OPS_H
