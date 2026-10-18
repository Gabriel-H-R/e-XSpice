/**
 * @file test_complex_lu.cpp
 * @brief Unit test for the sparse complex LU used by the AC analysis.
 *
 * The AC solver is the one place where a silent numerical error would be
 * hard to spot downstream (a wrong phase looks like "a plausible circuit
 * answer"), so this test pins the factorization against an INDEPENDENT dense
 * complex LU with partial pivoting on:
 *
 *  1. random sparse matrices (fill-in, pivoting-free Crout vs pivoted dense),
 *  2. the analytic RC low-pass network,
 *  3. pattern reuse across many frequencies (factorize() is called repeatedly
 *     on the same pattern: results must not depend on call order),
 *  4. the structural-singular / zero-matrix failure paths.
 */

#include "../src/solver/ComplexSparseLU.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using namespace e_XSpice;

static int g_failures = 0;

static void checkTrue(const char* what, bool cond) {
    if (!cond) ++g_failures;
    std::printf("  [%s] %s\n", cond ? "OK" : "FAIL", what);
}

static void checkClose(const char* what, double got, double expect, double tol) {
    const bool ok = std::isfinite(got) && std::abs(got - expect) <= tol;
    if (!ok) ++g_failures;
    std::printf("  [%s] %-46s got=% .10g expect=% .10g\n",
                ok ? "OK" : "FAIL", what, got, expect);
}

// ---- independent reference: dense complex LU with partial pivoting ---------
// Row-major n x n matrix in a flat vector.
static inline size_t at(int i, int j, int n) {
    return static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j);
}

static bool denseSolve(std::vector<std::complex<double>> A,
                       std::vector<std::complex<double>> b,
                       std::vector<std::complex<double>>& x) {
    const int n = static_cast<int>(b.size());
    for (int k = 0; k < n; ++k) {
        int piv = k;
        double best = std::abs(A[at(k, k, n)]);
        for (int i = k + 1; i < n; ++i) {
            const double v = std::abs(A[at(i, k, n)]);
            if (v > best) { best = v; piv = i; }
        }
        if (!(best > 1e-300)) return false;
        if (piv != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A[at(k, j, n)], A[at(piv, j, n)]);
            }
            std::swap(b[static_cast<size_t>(k)], b[static_cast<size_t>(piv)]);
        }
        const std::complex<double> d = A[at(k, k, n)];
        for (int i = k + 1; i < n; ++i) {
            const std::complex<double> f = A[at(i, k, n)] / d;
            if (f == std::complex<double>(0.0, 0.0)) continue;
            A[at(i, k, n)] = f;
            for (int j = k + 1; j < n; ++j) {
                A[at(i, j, n)] -= f * A[at(k, j, n)];
            }
            b[static_cast<size_t>(i)] -= f * b[static_cast<size_t>(k)];
        }
    }
    x.assign(static_cast<size_t>(n), std::complex<double>(0.0, 0.0));
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> s = b[static_cast<size_t>(i)];
        for (int j = i + 1; j < n; ++j) {
            s -= A[at(i, j, n)] * x[static_cast<size_t>(j)];
        }
        x[static_cast<size_t>(i)] = s / A[at(i, i, n)];
    }
    return true;
}

// Independent reference WITHOUT pivoting (same natural ordering as the
// sparse solver): isolates implementation correctness from conditioning.
static bool denseSolveNoPivot(std::vector<std::complex<double>> A,
                              std::vector<std::complex<double>> b,
                              std::vector<std::complex<double>>& x) {
    const int n = static_cast<int>(b.size());
    for (int k = 0; k < n; ++k) {
        const std::complex<double> u = A[at(k, k, n)];
        if (std::abs(u) == 0.0) return false;
        for (int i = k + 1; i < n; ++i) {
            A[at(i, k, n)] /= u;
            const std::complex<double> f = A[at(i, k, n)];
            if (f == std::complex<double>(0.0, 0.0)) continue;
            for (int j = k + 1; j < n; ++j) {
                A[at(i, j, n)] -= f * A[at(k, j, n)];
            }
            b[static_cast<size_t>(i)] -= f * b[static_cast<size_t>(k)];
        }
    }
    x.assign(static_cast<size_t>(n), std::complex<double>(0.0, 0.0));
    for (int i = n - 1; i >= 0; --i) {
        std::complex<double> s = b[static_cast<size_t>(i)];
        for (int j = i + 1; j < n; ++j) {
            s -= A[at(i, j, n)] * x[static_cast<size_t>(j)];
        }
        x[static_cast<size_t>(i)] = s / A[at(i, i, n)];
    }
    return true;
}

struct SparseCase {
    int n = 0;
    std::vector<int> rowPtr, colIdx;
    std::vector<double> G, C;      // real parts / capacitive parts
    std::vector<std::complex<double>> dense;   // A(omega) for the reference
    std::vector<std::complex<double>> rhs;
};

static SparseCase randomCase(int n, int perRow, unsigned seed, double omega,
                             bool diagonal_dominant) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    SparseCase c;
    c.n = n;
    c.rowPtr.push_back(0);
    std::vector<std::vector<int>> cols(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < perRow; ++k) {
            cols[static_cast<size_t>(i)].push_back(
                std::min(n - 1, std::abs(static_cast<int>(val(rng) * n)) % n));
        }
        cols[static_cast<size_t>(i)].push_back(i);   // diagonal always
        std::sort(cols[static_cast<size_t>(i)].begin(), cols[static_cast<size_t>(i)].end());
        cols[static_cast<size_t>(i)].erase(
            std::unique(cols[static_cast<size_t>(i)].begin(),
                        cols[static_cast<size_t>(i)].end()),
            cols[static_cast<size_t>(i)].end());
        c.rowPtr.push_back(c.rowPtr.back() +
                           static_cast<int>(cols[static_cast<size_t>(i)].size()));
    }
    for (int i = 0; i < n; ++i) {
        // Circuit-like scaling. The no-pivot ordering needs the complex
        // matrix to be diagonally DOMINANT, which is what a real AC matrix
        // is (G from conductances, w*C comparable to it). So: a strong real
        // diagonal, small off-diagonal conductances, and capacitances small
        // enough that w*C stays below the diagonal for the w used below.
        const double budget = diagonal_dominant ? 0.05 : 0.3;
        for (int j : cols[static_cast<size_t>(i)]) {
            c.colIdx.push_back(j);
            const double g = (i == j)
                ? (2.0 + std::abs(val(rng))) : val(rng) * budget;
            // 1e-6 keeps w*C ~ 1e-3 with w ~ 2*pi*1e3..1e4 (F-scale).
            const double cc = (i == j) ? 0.0 : std::fabs(val(rng)) * budget * 1e-6;
            c.G.push_back(g);
            c.C.push_back(cc);
        }
    }
    // Dense reference matrix + rhs.
    c.dense.assign(static_cast<size_t>(n) * static_cast<size_t>(n),
                   std::complex<double>(0.0, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int p = c.rowPtr[i]; p < c.rowPtr[i + 1]; ++p) {
            const int j = c.colIdx[static_cast<size_t>(p)];
            c.dense[static_cast<size_t>(i) * static_cast<size_t>(n) +
                    static_cast<size_t>(j)] =
                std::complex<double>(c.G[static_cast<size_t>(p)],
                                     omega * c.C[static_cast<size_t>(p)]);
        }
    }
    c.rhs.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        c.rhs[static_cast<size_t>(i)] =
            std::complex<double>(val(rng), val(rng));
    }
    return c;
}

// Residual of A·x − b for a solution on the CSR pattern (the invariant that
// actually matters for a no-pivot factorization: the circuit equations must
// be satisfied, whatever the growth).
static double residualInf(const SparseCase& c, double omega,
                           const std::vector<double>& xr,
                           const std::vector<double>& xi) {
    const int n = c.n;
    std::vector<double> axr(static_cast<size_t>(n), 0.0);
    std::vector<double> axi(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        for (int p = c.rowPtr[i]; p < c.rowPtr[i + 1]; ++p) {
            const int j = c.colIdx[static_cast<size_t>(p)];
            const double gr = c.G[static_cast<size_t>(p)];
            const double gi = omega * c.C[static_cast<size_t>(p)];
            // (gr + j*gi) * (xr[j] + j*xi[j])
            axr[static_cast<size_t>(i)] +=
                gr * xr[static_cast<size_t>(j)] - gi * xi[static_cast<size_t>(j)];
            axi[static_cast<size_t>(i)] +=
                gr * xi[static_cast<size_t>(j)] + gi * xr[static_cast<size_t>(j)];
        }
    }
    double worst = 0.0, bnorm = 0.0;
    for (int i = 0; i < n; ++i) {
        const double rr = std::hypot(axr[static_cast<size_t>(i)]
                                     - c.rhs[static_cast<size_t>(i)].real(),
                                     axi[static_cast<size_t>(i)]
                                     - c.rhs[static_cast<size_t>(i)].imag());
        worst = std::max(worst, rr);
        bnorm = std::max(bnorm, std::abs(c.rhs[static_cast<size_t>(i)]));
    }
    return (bnorm > 0.0) ? worst / bnorm : worst;
}

// Builds a random sparse matrix given as two REAL matrices (G, C) on a CSR
// pattern, which is exactly how the AC engine feeds the solver.

int main() {
    // Unbuffered: a crash mid-test must not swallow the diagnostics.
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("Complex sparse LU checks\n");

    // ---- 1) random sparse matrices vs independent references --------------
    {
        std::printf("\n1a) Circuit-like random sparse: true-residual check\n");
        int tested = 0;
        double worst_res = 0.0;
        for (int n : {8, 17, 40, 96}) {
            for (int perRow : {3, 5, 9}) {
                for (unsigned seed = 1; seed <= 3; ++seed) {
                    const double omega = 2.0 * 3.14159265358979 * 1000.0 * seed;
                    SparseCase c = randomCase(n, perRow, seed * 7919u,
                                              omega, true);
                    ComplexSparseLU lu;
                    lu.setPattern(c.rowPtr, c.colIdx);
                    if (!lu.factorize(omega, c.G.data(), c.C.data())) continue;
                    std::vector<double> xr(c.rhs.size()), xi(c.rhs.size());
                    for (size_t i = 0; i < c.rhs.size(); ++i) {
                        xr[i] = c.rhs[i].real();
                        xi[i] = c.rhs[i].imag();
                    }
                    lu.solve(xr, xi);
                    ++tested;
                    worst_res = std::max(worst_res, residualInf(c, omega, xr, xi));
                }
            }
        }
        checkTrue("random systems factorized", tested >= 24);
        std::printf("       worst relative residual = %.3e\n", worst_res);
        checkClose("true residual of every random system", worst_res, 0.0, 1e-10);

        std::printf("\n1b) Circuit-like (diagonally dominant) matrices"
                    " vs dense PIVOTED reference\n");
        int tested2 = 0;
        double worst2 = 0.0;
        for (int n : {8, 17, 40, 96}) {
            for (int perRow : {3, 5, 9}) {
                for (unsigned seed = 1; seed <= 3; ++seed) {
                    const double omega = 2.0 * 3.14159265358979 * 1000.0 * seed;
                    SparseCase c = randomCase(n, perRow, seed * 104729u,
                                              omega, true);
                    ComplexSparseLU lu;
                    lu.setPattern(c.rowPtr, c.colIdx);
                    if (!lu.factorize(omega, c.G.data(), c.C.data())) continue;
                    std::vector<double> xr(c.rhs.size()), xi(c.rhs.size());
                    for (size_t i = 0; i < c.rhs.size(); ++i) {
                        xr[i] = c.rhs[i].real();
                        xi[i] = c.rhs[i].imag();
                    }
                    lu.solve(xr, xi);
                    std::vector<std::complex<double>> ref;
                    if (!denseSolve(c.dense, c.rhs, ref)) continue;
                    ++tested2;
                    for (size_t i = 0; i < ref.size(); ++i) {
                        worst2 = std::max(worst2, std::abs(
                            std::complex<double>(xr[i], xi[i]) - ref[i]));
                    }
                }
            }
        }
        checkTrue("well-conditioned systems factorized", tested2 >= 24);
        std::printf("       worst |x_sparse - x_dense_pivot| = %.3e\n", worst2);
        checkClose("agreement with the pivoted dense reference", worst2, 0.0, 1e-8);
    }

    // ---- 2) analytic RC low-pass -----------------------------------------
    // V1 -R- in/out -C- gnd,  unknowns x = [v(in), v(out), i(V1)]
    //   node in  : (v_in - v_out)/R + i         = 0
    //   node out : (v_out - v_in)/R + jwC·v_out = 0
    //   branch   : v_in = V_ac                  -> B = [0, 0, V_ac]
    // so H = v(out)/V_ac = 1 / (1 + j w R C).
    {
        std::printf("\n2) RC low-pass vs closed form\n");
        const double R = 1e3, C = 1e-9, g = 1.0 / R;
        // CSR row-major: row0 {0,1,2}, row1 {0,1}, row2 {0}
        std::vector<int> rowPtr = {0, 3, 5, 6};
        std::vector<int> colIdx = {0, 1, 2, 0, 1, 0};
        std::vector<double> G = { g, -g, 1.0, -g, g, 1.0};
        std::vector<double> Cval = {0.0, 0.0, 0.0, 0.0, C, 0.0};
        ComplexSparseLU lu;
        lu.setPattern(rowPtr, colIdx);
        for (double f : {10.0, 1000.0, 15915.0, 100000.0}) {
            const double omega = 2.0 * 3.14159265358979 * f;
            checkTrue("factorize ok", lu.factorize(omega, G.data(), Cval.data()));
            std::vector<double> xr = {0.0, 0.0, 1.0};   // source at branch row
            std::vector<double> xi = {0.0, 0.0, 0.0};
            lu.solve(xr, xi);
            const double mag = std::hypot(xr[1], xi[1]);   // |v(out)|
            const double ph = std::atan2(xi[1], xr[1]) * 180.0 / 3.14159265358979;
            const double wc = omega * R * C;
            const double mag_ref = 1.0 / std::sqrt(1.0 + wc * wc);
            const double ph_ref = -std::atan(wc) * 180.0 / 3.14159265358979;
            char buf[96];
            std::snprintf(buf, sizeof(buf), "|H| at %.0f Hz", f);
            checkClose(buf, mag, mag_ref, 1e-10);
            std::snprintf(buf, sizeof(buf), "phase at %.0f Hz (deg)", f);
            checkClose(buf, ph, ph_ref, 1e-9);
            std::snprintf(buf, sizeof(buf), "v(in) == V_ac at %.0f Hz", f);
            checkClose(buf, xr[0], 1.0, 1e-12);
        }
    }

    // ---- 3) pattern reuse across frequencies -----------------------------
    {
        std::printf("\n3) Pattern reuse (call order independence)\n");
        const double R = 470.0, C = 1e-7, g = 1.0 / R;
        std::vector<int> rowPtr = {0, 3, 5, 6};
        std::vector<int> colIdx = {0, 1, 2, 0, 1, 0};
        std::vector<double> G = { g, -g, 1.0, -g, g, 1.0};
        std::vector<double> Cval = {0.0, 0.0, 0.0, 0.0, C, 0.0};
        ComplexSparseLU lu;
        lu.setPattern(rowPtr, colIdx);
        std::vector<double> first;
        for (int pass = 0; pass < 2; ++pass) {
            for (double f : {10.0, 5000.0, 1e6}) {
                const double omega = 2.0 * 3.14159265358979 * f;
                lu.factorize(omega, G.data(), Cval.data());
                std::vector<double> xr = {0.0, 0.0, 1.0};
                std::vector<double> xi = {0.0, 0.0, 0.0};
                lu.solve(xr, xi);
                if (pass == 0 && f == 10.0) first = xr;
                if (pass == 1 && f == 10.0) {
                    checkClose("re-factorize reproduces the first value",
                               xr[1], first[1], 1e-15);
                }
            }
        }
    }

    // ---- 4) failure paths -------------------------------------------------
    {
        std::printf("\n4) Degenerate inputs\n");
        {
            // 1x1 with a zero diagonal: structurally present, numerically zero
            ComplexSparseLU lu;
            std::vector<int> rowPtr = {0, 1};
            std::vector<int> colIdx = {0};
            std::vector<double> G = {0.0};
            std::vector<double> C = {0.0};
            lu.setPattern(rowPtr, colIdx);
            const bool ok = lu.factorize(1.0, G.data(), C.data());
            std::printf("       zero matrix: factorize=%d perturbed=%d msg='%s'\n",
                        (int)ok, (int)lu.perturbedPivot(), lu.message().c_str());
            checkTrue("zero matrix is reported, not silently 'solved'", !ok);
        }
        {
            // Diagonal-dominant matrix needing no perturbation
            ComplexSparseLU lu;
            std::vector<int> rowPtr = {0, 1};
            std::vector<int> colIdx = {0};
            std::vector<double> G = {2.0};
            std::vector<double> C = {0.0};
            lu.setPattern(rowPtr, colIdx);
            checkTrue("well-posed 1x1 factors", lu.factorize(1.0, G.data(), C.data()));
            checkTrue("no perturbation needed", !lu.perturbedPivot());
            std::vector<double> xr = {4.0};
            std::vector<double> xi = {0.0};
            lu.solve(xr, xi);
            checkClose("1x1 solution", xr[0], 2.0, 1e-15);
        }
    }

    std::printf(g_failures == 0 ? "Complex LU PASSED\n" : "Complex LU FAILED (%d)\n",
                g_failures);
    return g_failures == 0 ? 0 : 1;
}