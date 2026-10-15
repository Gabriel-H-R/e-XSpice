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

SolverResult LinearSolver::solve(const MatrixOps& A,
                                 const std::vector<double>& b,
                                 std::vector<double>& x,
                                 const std::vector<double>* x0) {
    if (A.rows() <= 0) {
        throw std::runtime_error("Solver received an empty matrix");
    }
    if (static_cast<int>(b.size()) != A.rows()) {
        throw std::runtime_error("Vector size must match matrix dimensions");
    }
    
    // Choose solver based on type
    switch (type_) {
        case SolverType::CONJUGATE_GRADIENT:
            return solveConjugateGradient(A, b, x, x0);
        case SolverType::BICGSTAB:
            return solveBiCGSTAB(A, b, x, x0);
        case SolverType::DIRECT_LU:
            return solveDenseLU(A, b, x);
        default:
            throw std::runtime_error("Solver type not implemented");
    }
}

double LinearSolver::dotProduct(const std::vector<double>& a,
                                const std::vector<double>& b) const {
    const int n = static_cast<int>(a.size());
    constexpr int kParallelMin = 512;   // below this, serial is faster
    if (n >= kParallelMin) {
        double result = 0.0;
        #pragma omp parallel for reduction(+:result)
        for (int i = 0; i < n; ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double LinearSolver::vectorNorm(const std::vector<double>& v) const {
    return std::sqrt(dotProduct(v, v));
}

bool LinearSolver::verifyTrueResidual(const MatrixOps& A,
                                      const std::vector<double>& b,
                                      const std::vector<double>& x,
                                      double b_norm,
                                      SolverResult& result) const {
    std::vector<double> Ax;
    A.multiplyVector(x, Ax);
    double sum = 0.0;
    const size_t n = b.size();
    for (size_t i = 0; i < n; ++i) {
        const double d = b[i] - Ax[i];
        sum += d * d;
    }
    const double true_norm = std::sqrt(sum);
    result.residual_norm = true_norm;
    if (true_norm / b_norm < tolerance_) {
        result.converged = true;
        result.message = "Converged";
        return true;
    }
    result.converged = false;
    result.message = "True residual too large (recurrence drifted)";
    return false;
}

void LinearSolver::vectorScale(std::vector<double>& v, double alpha) const {
    const int n = static_cast<int>(v.size());

    if (n >= 512) {
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            v[i] *= alpha;
        }
        return;
    }
    for (int i = 0; i < n; ++i) {
        v[i] *= alpha;
    }
}

void LinearSolver::vectorAdd(std::vector<double>& result,
                             const std::vector<double>& a,
                             double alpha,
                             const std::vector<double>& b) const {
    const int n = static_cast<int>(a.size());

    if (n >= 512) {
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            result[i] = a[i] + alpha * b[i];
        }
        return;
    }
    for (int i = 0; i < n; ++i) {
        result[i] = a[i] + alpha * b[i];
    }
}

SolverResult LinearSolver::solveConjugateGradient(
    const MatrixOps& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const std::vector<double>* x0) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int n = A.rows();
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
    const MatrixOps& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const std::vector<double>* x0) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    const int n = A.rows();
    
    // Initial guess
    x.assign(n, 0.0);
    if (x0 != nullptr && static_cast<int>(x0->size()) == n) {
        x = *x0;
    }
    
    // Jacobi (diagonal) preconditioner: cheap and very effective for MNA
    // systems, whose near-singular rows come from diode exponentials.
    std::vector<double> minv(n, 1.0);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        const double d = A.get(i, i);
        minv[i] = (std::abs(d) > 1e-300) ? 1.0 / d : 1.0;
    }
    
    // r = b − A·x
    std::vector<double> Ax(n);
    std::vector<double> r(n);
    A.multiplyVector(x, Ax);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];
    
    std::vector<double> rhat = r;            // shadow residual
    std::vector<double> p(n, 0.0), v(n, 0.0);
    std::vector<double> s(n), t(n), ph(n), sh(n);
    
    double rho_prev = 1.0, alpha = 1.0, omega = 1.0;
    double b_norm = vectorNorm(b);
    if (b_norm < 1e-300) b_norm = 1.0;
    
    SolverResult result;
    result.converged = false;
    result.iterations = 0;
    
    auto converged = [&](const std::vector<double>& resid) {
        const double rn = vectorNorm(resid);
        result.residual_norm = rn;
        return rn / b_norm < tolerance_;
    };
    
    if (converged(r)) {                      // already a solution
        result.converged = true;
        result.message = "Converged";
        result.iterations = 0;
        return result;
    }
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        result.iterations = iter + 1;
        
        // ρ = r̂·r, with a safeguard restart when it vanishes.
        double rho = dotProduct(rhat, r);
        int restarts = 0;
        while (std::abs(rho) < 1e-300 && restarts < 2) {
            rhat = r;                          // re-orthogonalize / restart
            rho = dotProduct(rhat, r);
            ++restarts;
        }
        if (std::abs(rho) < 1e-300) {
            result.message = "BiCGSTAB breakdown: rho near zero";
            break;
        }
        
        double beta = (iter == 0) ? 0.0
                                  : (rho / rho_prev) * (alpha / omega);
        // p = r + β·(p − ω·v)
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        
        // p̂ = M⁻¹·p ; v = A·p̂
        for (int i = 0; i < n; ++i) ph[i] = minv[i] * p[i];
        A.multiplyVector(ph, v);
        
        const double denom = dotProduct(rhat, v);
        if (std::abs(denom) < 1e-300) {
            result.message = "BiCGSTAB breakdown: r^T·A·p near zero";
            break;
        }
        alpha = rho / denom;
        
        // s = r − α·v  (early convergence test on ‖s‖)
        for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
        if (converged(s)) {
            for (int i = 0; i < n; ++i) x[i] += alpha * ph[i];
            // Guard against recurrence drift: only a true residual counts.
            // On failure the caller falls back to the direct solver.
            verifyTrueResidual(A, b, x, b_norm, result);
            break;
        }
        
        // ŝ = M⁻¹·s ; t = A·ŝ
        for (int i = 0; i < n; ++i) sh[i] = minv[i] * s[i];
        A.multiplyVector(sh, t);
        const double tt = dotProduct(t, t);
        if (std::abs(tt) < 1e-300) {
            result.message = "BiCGSTAB breakdown: t^T·t near zero";
            break;
        }
        omega = dotProduct(t, s) / tt;
        
        // x += α·p̂ + ω·ŝ
        for (int i = 0; i < n; ++i) x[i] += alpha * ph[i] + omega * sh[i];

        // r = s − ω·t
        for (int i = 0; i < n; ++i) r[i] = s[i] - omega * t[i];

        if (converged(r)) {
            // Guard against recurrence drift (see above).
            verifyTrueResidual(A, b, x, b_norm, result);
            break;
        }
        if (std::abs(omega) < 1e-300) {
            result.message = "BiCGSTAB breakdown: omega near zero";
            break;
        }
        rho_prev = rho;
        
        if (verbose_ && iter % 100 == 0) {
            std::cout << "  BiCGSTAB iteration " << iter
                      << ", residual = " << result.residual_norm << std::endl;
        }
    }
    
    if (!result.converged && result.message.empty()) {
        result.message = "Did not converge within max iterations";
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.elapsed_time =
        std::chrono::duration<double>(end_time - start_time).count();
    return result;
}

SolverResult LinearSolver::solveDenseLU(const MatrixOps& A,
                                        const std::vector<double>& b,
                                        std::vector<double>& x) const {
    auto start_time = std::chrono::high_resolution_clock::now();

    const int n = A.rows();
    x.assign(n, 0.0);

    // Gather into a dense row-major workspace (LU is done in place).
    std::vector<double> lu(static_cast<size_t>(n) * n, 0.0);
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n; ++i) {
        double* row = &lu[static_cast<size_t>(i) * n];
        for (int j = 0; j < n; ++j) {
            row[j] = A.get(i, j);
        }
    }

    SolverResult result;
    result.converged = false;

    // Doolittle LU with partial pivoting on rows.
    std::vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;

    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double best = std::abs(lu[static_cast<size_t>(k) * n + k]);
        for (int i = k + 1; i < n; ++i) {
            const double v = std::abs(lu[static_cast<size_t>(i) * n + k]);
            if (v > best) { best = v; pivot = i; }
        }
        if (best < 1e-300) {
            result.message = "LU: matrix is singular";
            result.iterations = 0;
            result.elapsed_time =
                std::chrono::duration<double>(
                    std::chrono::high_resolution_clock::now() - start_time).count();
            return result;
        }
        if (pivot != k) {
            #pragma omp parallel for
            for (int j = 0; j < n; ++j) {
                std::swap(lu[static_cast<size_t>(k) * n + j],
                          lu[static_cast<size_t>(pivot) * n + j]);
            }
            std::swap(perm[k], perm[pivot]);
        }

        const double akk = lu[static_cast<size_t>(k) * n + k];
        #pragma omp parallel for schedule(dynamic, 8)
        for (int i = k + 1; i < n; ++i) {
            const double factor = lu[static_cast<size_t>(i) * n + k] / akk;
            lu[static_cast<size_t>(i) * n + k] = factor;
            if (factor == 0.0) continue;
            double* row_i = &lu[static_cast<size_t>(i) * n];
            const double* row_k = &lu[static_cast<size_t>(k) * n];
            for (int j = k + 1; j < n; ++j) {
                row_i[j] -= factor * row_k[j];
            }
        }
    }

    // Forward substitution (L·y = P·b), then back substitution (U·x = y).
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = b[perm[i]];
        const double* row = &lu[static_cast<size_t>(i) * n];
        for (int j = 0; j < i; ++j) sum -= row[j] * y[j];
        y[i] = sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[i];
        const double* row = &lu[static_cast<size_t>(i) * n];
        for (int j = i + 1; j < n; ++j) sum -= row[j] * x[j];
        x[i] = sum / row[i];
    }

    // Residual check.
    std::vector<double> Ax;
    A.multiplyVector(x, Ax);
    double resid = 0.0, bnorm = 0.0;
    for (int i = 0; i < n; ++i) {
        resid += (b[i] - Ax[i]) * (b[i] - Ax[i]);
        bnorm += b[i] * b[i];
    }
    resid = std::sqrt(resid);
    bnorm = std::sqrt(bnorm);
    result.residual_norm = resid;
    result.iterations = 1;
    result.converged = (resid < 1e-8 * std::max(1.0, bnorm));
    result.message = result.converged ? "Converged" : "LU residual too large";
    result.elapsed_time =
        std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - start_time).count();
    return result;
}

} // namespace e_XSpice
