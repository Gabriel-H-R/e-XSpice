/**
 * @author Gabriel Henrique Silva
 * @file LinearSolver.cpp
 * @brief Implementation of iterative linear solvers
 */

#include "LinearSolver.h"
#include <iostream>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <omp.h>

namespace e_XSpice {

LinearSolver::LinearSolver(SolverType type)
    : type_(type)
    , max_iterations_(10000)
    , tolerance_(1e-9)
    , verbose_(false)
{
}

SolverResult LinearSolver::solve(const SparseMatrix& A,
                                 const std::vector<double>& b,
                                 std::vector<double>& x,
                                 const std::vector<double>* x0) {
    if (A.getFormat() != MatrixFormat::CSR) {
        throw std::runtime_error("Matrix must be in CSR format");
    }
    
    if (A.getRows() != A.getCols()) {
        throw std::runtime_error("Matrix must be square");
    }
    
    if (static_cast<int>(b.size()) != A.getRows()) {
        throw std::runtime_error("Vector size must match matrix dimensions");
    }
    
    // Choose solver based on type
    switch (type_) {
        case SolverType::CONJUGATE_GRADIENT:
            return solveConjugateGradient(A, b, x, x0);
        case SolverType::BICGSTAB:
            return solveBiCGSTAB(A, b, x, x0);
        default:
            throw std::runtime_error("Solver type not implemented");
    }
}

double LinearSolver::dotProduct(const std::vector<double>& a,
                                const std::vector<double>& b) const {
    int n = a.size();
    double result = 0.0;
    
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < n; ++i) {
        result += a[i] * b[i];
    }
    
    return result;
}

double LinearSolver::vectorNorm(const std::vector<double>& v) const {
    return std::sqrt(dotProduct(v, v));
}

void LinearSolver::vectorScale(std::vector<double>& v, double alpha) const {
    int n = v.size();
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        v[i] *= alpha;
    }
}

void LinearSolver::vectorAdd(std::vector<double>& result,
                             const std::vector<double>& a,
                             double alpha,
                             const std::vector<double>& b) const {
    int n = a.size();
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        result[i] = a[i] + alpha * b[i];
    }
}

SolverResult LinearSolver::solveConjugateGradient(
    const SparseMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const std::vector<double>* x0) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int n = A.getRows();
    x.resize(n, 0.0);
    
    // Initialization
    if (x0 != nullptr) {
        x = *x0;
    }
    
    // r = b - Ax
    std::vector<double> r(n);
    std::vector<double> Ax(n);
    A.multiplyVector(x, Ax);
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax[i];
    }
    
    std::vector<double> p = r;  // Initial search direction
    double rsold = dotProduct(r, r);
    
    SolverResult result;
    result.converged = false;
    result.iterations = 0;
    
    double b_norm = vectorNorm(b);
    if (b_norm < tolerance_) {
        b_norm = 1.0;  // Avoids division by zero
    }
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        result.iterations = iter + 1;
        
        // Ap = A * p
        std::vector<double> Ap(n);
        A.multiplyVector(p, Ap);
        
        // alpha = rsold / (p' * Ap)
        double pAp = dotProduct(p, Ap);
        if (std::abs(pAp) < 1e-30) {
            result.message = "Algorithm breakdown: p'Ap near zero";
            break;
        }
        double alpha = rsold / pAp;
        
        // x = x + alpha * p
        vectorAdd(x, x, alpha, p);
        
        // r = r - alpha * Ap
        vectorAdd(r, r, -alpha, Ap);
        
        double rsnew = dotProduct(r, r);
        result.residual_norm = std::sqrt(rsnew);
        
        // Convergence test (relative residual)
        if (result.residual_norm / b_norm < tolerance_) {
            result.converged = true;
            result.message = "Converged";
            break;
        }
        
        // beta = rsnew / rsold
        double beta = rsnew / rsold;
        
        // p = r + beta * p
        vectorAdd(p, r, beta, p);
        
        rsold = rsnew;
        
        if (verbose_ && (iter % 100 == 0)) {
            std::cout << "  CG Iteration " << iter 
                     << ", residual = " << result.residual_norm << std::endl;
        }
    }
    
    if (!result.converged) {
        result.message = "Did not converge within max iterations";
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.elapsed_time = std::chrono::duration<double>(end_time - start_time).count();
    
    return result;
}

SolverResult LinearSolver::solveBiCGSTAB(
    const SparseMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const std::vector<double>* x0) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int n = A.getRows();
    x.resize(n, 0.0);
    
    // Initialization
    if (x0 != nullptr) {
        x = *x0;
    }
    
    // r = b - Ax
    std::vector<double> r(n);
    std::vector<double> Ax(n);
    A.multiplyVector(x, Ax);
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax[i];
    }
    
    std::vector<double> r0 = r;  // r0 arbitrary (can be initial r)
    std::vector<double> p = r;
    std::vector<double> v(n, 0.0);
    
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    
    SolverResult result;
    result.converged = false;
    result.iterations = 0;
    
    double b_norm = vectorNorm(b);
    if (b_norm < tolerance_) {
        b_norm = 1.0;
    }
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        result.iterations = iter + 1;
        
        double rho_new = dotProduct(r0, r);
        
        if (std::abs(rho_new) < 1e-30) {
            result.message = "BiCGSTAB breakdown: rho near zero";
            break;
        }
        
        if (iter > 0) {
            double beta = (rho_new / rho) * (alpha / omega);
            
            // p = r + beta * (p - omega * v)
            std::vector<double> temp(n);
            vectorAdd(temp, p, -omega, v);
            vectorAdd(p, r, beta, temp);
        }
        
        // v = A * p
        A.multiplyVector(p, v);
        
        double r0v = dotProduct(r0, v);
        if (std::abs(r0v) < 1e-30) {
            result.message = "BiCGSTAB breakdown: r0'v near zero";
            break;
        }
        
        alpha = rho_new / r0v;
        
        // s = r - alpha * v
        std::vector<double> s(n);
        vectorAdd(s, r, -alpha, v);
        
        // Early convergence test
        result.residual_norm = vectorNorm(s);
        if (result.residual_norm / b_norm < tolerance_) {
            vectorAdd(x, x, alpha, p);
            result.converged = true;
            result.message = "Converged";
            break;
        }
        
        // t = A * s
        std::vector<double> t(n);
        A.multiplyVector(s, t);
        
        double tt = dotProduct(t, t);
        if (tt < 1e-30) {
            result.message = "BiCGSTAB breakdown: t't near zero";
            break;
        }
        
        omega = dotProduct(t, s) / tt;
        
        // x = x + alpha * p + omega * s
        vectorAdd(x, x, alpha, p);
        vectorAdd(x, x, omega, s);
        
        // r = s - omega * t
        vectorAdd(r, s, -omega, t);
        
        result.residual_norm = vectorNorm(r);
        
        if (result.residual_norm / b_norm < tolerance_) {
            result.converged = true;
            result.message = "Converged";
            break;
        }
        
        if (std::abs(omega) < 1e-30) {
            result.message = "BiCGSTAB breakdown: omega near zero";
            break;
        }
        
        rho = rho_new;
        
        if (verbose_ && (iter % 100 == 0)) {
            std::cout << "  BiCGSTAB Iteration " << iter 
                     << ", residual = " << result.residual_norm << std::endl;
        }
    }
    
    if (!result.converged) {
        result.message = "Did not converge within max iterations";
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.elapsed_time = std::chrono::duration<double>(end_time - start_time).count();
    
    return result;
}

} // namespace e_XSpice
