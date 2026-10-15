/**
 * @author Gabriel Henrique Silva
 * @file DistributedKrylov.cpp
 * @brief Row-distributed Jacobi-preconditioned BiCGSTAB (MPI).
 */

#include "DistributedKrylov.h"
#include "MPIManager.h"
#include "../solver/DaeMatrix.h"
#include "../solver/LinearSolver.h"

#include <omp.h>
#include <cmath>
#include <chrono>
#include <algorithm>

namespace e_XSpice {

namespace {

// Partial dot over a row block (serial when tiny, matching solver policy).
double blockDot(const std::vector<double>& a, const std::vector<double>& b,
                int r0, int r1) {
    const int count = r1 - r0;
    if (count <= 0) return 0.0;
    double sum = 0.0;
    if (count >= 512) {
        #pragma omp parallel for reduction(+:sum)
        for (int i = r0; i < r1; ++i) sum += a[i] * b[i];
    } else {
        for (int i = r0; i < r1; ++i) sum += a[i] * b[i];
    }
    return sum;
}

double fullNorm(const std::vector<double>& v, int r0, int r1, MPIManager& mpi) {
    return std::sqrt(mpi.allReduceSumDouble(blockDot(v, v, r0, r1)));
}

// True residual from scratch, computed REDUNDANTLY (A, x, b are identical on
// every rank): zero communication, identical verdict everywhere. Guards the
// same recurrence-drift failure as the serial path (see LinearSolver).
bool verifyTrueResidualDist(const DaeMatrix& A,
                            const std::vector<double>& b,
                            const std::vector<double>& x,
                            double b_norm, double tolerance,
                            SolverResult& result) {
    std::vector<double> Ax;
    A.multiplyVector(x, Ax);
    double sum = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        const double d = b[i] - Ax[i];
        sum += d * d;
    }
    const double true_norm = std::sqrt(sum);
    result.residual_norm = true_norm;
    if (true_norm / b_norm < tolerance) {
        result.converged = true;
        result.message = "Converged";
        return true;
    }
    result.converged = false;
    result.message = "True residual too large (recurrence drifted)";
    return false;
}

} // namespace

SolverResult distributedBiCGSTAB(const DaeMatrix& A,
                                 const std::vector<double>& b,
                                 std::vector<double>& x,
                                 int rowBegin, int rowEnd,
                                 MPIManager& mpi,
                                 int maxIterations,
                                 double tolerance) {
    auto start_time = std::chrono::high_resolution_clock::now();
    SolverResult result;
    result.converged = false;
    result.iterations = 0;

    const int n = A.rows();
    if (static_cast<int>(b.size()) != n) {
        result.message = "distributed BiCGSTAB: size mismatch";
        return result;
    }
    if (static_cast<int>(x.size()) != n) {
        x.assign(n, 0.0);
    }

    // Jacobi preconditioner (replicated, identical everywhere).
    std::vector<double> minv(n, 1.0);
    for (int i = 0; i < n; ++i) {
        const double d = A.get(i, i);
        minv[i] = (std::abs(d) > 1e-30) ? 1.0 / d : 1.0;
    }

    std::vector<double> r(n), rhat(n), p(n, 0.0), v(n, 0.0);
    std::vector<double> s(n), t(n), ph(n), sh(n), tmp(n, 0.0);

    // r = b − A·x (block matvec + fuse; x is the replicated guess).
    A.multiplyVectorBlock(x, tmp, rowBegin, rowEnd);
    std::vector<double> Ax(n, 0.0);
    mpi.allGatherV(
        std::vector<double>(tmp.begin() + rowBegin, tmp.begin() + rowEnd),
        Ax, n);
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ax[i];
    rhat = r;

    double b_norm = fullNorm(b, rowBegin, rowEnd, mpi);
    if (b_norm < 1e-300) b_norm = 1.0;

    if (fullNorm(r, rowBegin, rowEnd, mpi) / b_norm < tolerance) {
        result.converged = true;
        result.message = "Converged";
        result.residual_norm = 0.0;
        result.elapsed_time =
            std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - start_time).count();
        return result;
    }

    double rho_prev = 1.0, alpha = 1.0, omega = 1.0;

    for (int iter = 0; iter < maxIterations; ++iter) {
        result.iterations = iter + 1;

        // Breakdown guards mirror the serial BiCGSTAB exactly (1e-30), so
        // both paths fall back to the direct solver on the same systems and
        // Newton sees identical behavior. (A looser guard lets a
        // near-zero denominator through, producing a huge alpha/omega and a
        // spurious "converged" correction — observed on ~1e-18 RHS norms.)
        double rho = mpi.allReduceSumDouble(blockDot(rhat, r, rowBegin, rowEnd));
        if (std::abs(rho) < 1e-30) {
            result.message = "distributed BiCGSTAB breakdown: rho near zero";
            break;
        }

        const double beta = (iter == 0) ? 0.0 : (rho / rho_prev) * (alpha / omega);
        for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);

        for (int i = 0; i < n; ++i) ph[i] = minv[i] * p[i];
        A.multiplyVectorBlock(ph, tmp, rowBegin, rowEnd);
        mpi.allGatherV(
            std::vector<double>(tmp.begin() + rowBegin, tmp.begin() + rowEnd),
            v, n);

        const double denom =
            mpi.allReduceSumDouble(blockDot(rhat, v, rowBegin, rowEnd));
        if (std::abs(denom) < 1e-30) {
            result.message = "distributed BiCGSTAB breakdown: r^T·A·p near zero";
            break;
        }
        alpha = rho / denom;

        for (int i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
        if (fullNorm(s, rowBegin, rowEnd, mpi) / b_norm < tolerance) {
            for (int i = 0; i < n; ++i) x[i] += alpha * ph[i];
            // Guard against recurrence drift (falls back to direct solve).
            verifyTrueResidualDist(A, b, x, b_norm, tolerance, result);
            break;
        }

        for (int i = 0; i < n; ++i) sh[i] = minv[i] * s[i];
        A.multiplyVectorBlock(sh, tmp, rowBegin, rowEnd);
        mpi.allGatherV(
            std::vector<double>(tmp.begin() + rowBegin, tmp.begin() + rowEnd),
            t, n);

        // Batch the two remaining dots (t·t and t·s) into ONE Allreduce of
        // 2 doubles — collectives are latency-bound, so halving the round
        // count matters more than halving the payload.
        const std::vector<double> dots_in = {
            blockDot(t, t, rowBegin, rowEnd),
            blockDot(t, s, rowBegin, rowEnd)};
        std::vector<double> dots_out;
        mpi.allReduceSum(dots_in, dots_out);
        const double tt = dots_out[0];
        const double ts = dots_out[1];
        if (std::abs(tt) < 1e-30) {
            result.message = "distributed BiCGSTAB breakdown: t^T·t near zero";
            break;
        }
        omega = ts / tt;

        for (int i = 0; i < n; ++i) x[i] += alpha * ph[i] + omega * sh[i];
        for (int i = 0; i < n; ++i) r[i] = s[i] - omega * t[i];

        if (fullNorm(r, rowBegin, rowEnd, mpi) / b_norm < tolerance) {
            // Guard against recurrence drift (falls back to direct solve).
            verifyTrueResidualDist(A, b, x, b_norm, tolerance, result);
            break;
        }
        if (std::abs(omega) < 1e-30) {
            result.message = "distributed BiCGSTAB breakdown: omega near zero";
            break;
        }
        rho_prev = rho;
    }

    if (!result.converged && result.message.empty()) {
        result.message = "Did not converge within max iterations";
    }
    result.residual_norm = fullNorm(r, rowBegin, rowEnd, mpi);
    result.elapsed_time =
        std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - start_time).count();
    return result;
}

} // namespace e_XSpice
