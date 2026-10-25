/**
 * @author Gabriel Henrique Silva
 * @file HBAnalysis.cpp
 * @brief Harmonic Balance analysis (`.HB`).
 *
 * See HBAnalysis.h for the formulation. The engine:
 *
 *   1. runs the operating point (the DC harmonic) and, when the netlist
 *      carries an `AC` excitation, the AC solution at the fundamental (the
 *      first harmonic) — the initial guess, without a transient run;
 *   2. iterates Newton on the harmonic system:
 *        a. reconstructs x(t_n) at the phase points of one period;
 *        b. evaluates the devices there (F, Q, B) and their Jacobians;
 *        c. transforms the residual and the Jacobian back to the frequency
 *           domain (the DFT);
 *        d. solves the expanded harmonic system and updates the harmonics;
 *      with a step damping (a halved step) when the plain Newton step does
 *      not reduce the residual — the harmonic balance Newton is famously
 *      prone to overshooting on a strongly nonlinear circuit.
 *   3. reports the harmonic coefficients and the reconstructed waveform.
 *
 * INDEXING (the part that is easy to get wrong)
 * ----------------------------------------------
 *   - the harmonics of the SOLUTION are k = 0 … K (k = 0 is the DC);
 *   - the expanded system carries the blocks b = 0 … 2K, the block b being
 *     the harmonic p_b = b − K (so the negative harmonics are the
 *     conjugates of the positive ones);
 *   - the DFT of the time-domain Jacobian is needed for the coefficients
 *     p = −2K … 2K (the intermodulation products), so the DFT slots are
 *     q = 0 … 4K with p = q − 2K, and the block (b, m) reads the slot
 *     q = (b − m) + 2K;
 *   - the number of phase points must exceed 4K, or the highest
 *     intermodulation coefficients alias onto the others.
 */

#include "HBAnalysis.h"

#include "../core/DaeSolution.h"
#include "../core/Device.h"
#include "../core/Node.h"
#include "../solver/ComplexSparseLU.h"
#include "../solver/DaeMatrix.h"
#include "DaeSystem.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace e_XSpice {

namespace {

constexpr double kTwoPi = 6.28318530717958647692;

} // namespace

HBAnalysis::HBAnalysis(Circuit& circuit) : circuit_(circuit) {}
HBAnalysis::~HBAnalysis() = default;

HBResult HBAnalysis::run(bool use_mpi) {
    HBResult result;
    result.fundamental = fundamental_;
    result.num_harmonics = harmonics_;

    // The harmonic solve is a dense block problem: the operating point and
    // the device evaluation are replicated on every rank, so the MPI layout
    // would only add collectives without changing the result.
    (void)use_mpi;

    if (fundamental_ <= 0.0) {
        result.message = "the fundamental frequency must be positive";
        return result;
    }
    if (harmonics_ < 1) {
        result.message = "at least one harmonic is required";
        return result;
    }

    const int K = harmonics_;
    const int nb = 2 * K + 1;      // the blocks of the expanded system
    const int nq = 4 * K + 1;      // the DFT slots (the harmonics -2K … 2K)

    // ---- the DAE layout ----------------------------------------------------
    DaeSystem system(circuit_);
    system.setup();
    const int n = system.size();
    const auto& node_map = system.nodeMap();
    result.num_unknowns = n;
    result.node_map = node_map;

    // The HB evaluates the devices at many points: the Newton limiting would
    // change the residual it is trying to solve, so it is turned off.
    for (const auto& d : circuit_.getDevices()) {
        d->setLimitingEnabled(false);
    }

    // ---- the operating point (the DC harmonic) -----------------------------
    system.setMaxIterations(max_iterations_);
    system.setRelTol(1e-3);
    system.setAbsTol(1e-12);
    system.setGminSchedule(1e-3, 1e-12);
    system.setSolverType(SolverType::BICGSTAB);
    system.beginSolve(/*initial_solve=*/true);
    {
        const NLSolveResult op = system.solveNonlinear(/*time=*/0.0,
                                                      /*q_coeff=*/0.0,
                                                      /*use_history=*/false);
        if (!op.converged) {
            result.message = "the operating point did not converge: " + op.message;
            return result;
        }
        system.extractSolution();
    }

    // ---- the unknowns: the harmonic coefficients X_k ----------------------
    std::vector<std::complex<double>> X(static_cast<size_t>(K + 1) *
                                         static_cast<size_t>(n), {0.0, 0.0});
    for (int i = 0; i < n; ++i) {
        X[static_cast<size_t>(i)] = system.solution()[static_cast<size_t>(i)];
    }

    // The first harmonic: the AC solution at the fundamental (linearized
    // about the operating point). It is a good initial guess and costs one
    // factorization — much cheaper than Xyce's transient-assisted start.
    {
        const double omega = kTwoPi * fundamental_;
        system.refreshLinearization(/*include_gmin=*/false);
        const DaeMatrix& G = system.jacobianF();
        const DaeMatrix& C = system.jacobianQ();
        ComplexSparseLU lu;
        lu.setPattern(G.rowPointers(), G.colIndices());
        if (lu.factorize(omega, G.values().data(), C.values().data())) {
            std::vector<double> bre(static_cast<size_t>(n), 0.0);
            std::vector<double> bim(static_cast<size_t>(n), 0.0);
            for (const auto& d : circuit_.getDevices()) {
                d->loadACExcitation(omega, bre, bim);
            }
            lu.solve(bre, bim);
            for (int i = 0; i < n; ++i) {
                X[static_cast<size_t>(n) + static_cast<size_t>(i)] =
                    std::complex<double>(bre[static_cast<size_t>(i)],
                                        bim[static_cast<size_t>(i)]);
            }
        }
    }

    // ---- the time grid -----------------------------------------------------
    // N samples per period of the fundamental; N must exceed 4K so the
    // intermodulation coefficients do not alias onto each other.
    int N = max_time_points_ > 0 ? max_time_points_ : (4 * K + 1);
    if (N < 4 * K + 1) N = 4 * K + 1;
    result.num_time_points = N;
    const double period = 1.0 / fundamental_;
    const double domega = kTwoPi / static_cast<double>(N);
    const double omega0 = kTwoPi * fundamental_;

    // The DFT weights e^{-j·2π·t·q/N} and the synthesis weights
    // e^{+j·2π·t·m/N}.
    std::vector<std::complex<double>> w_dft(
        static_cast<size_t>(N) * static_cast<size_t>(nq), {0.0, 0.0});
    std::vector<std::complex<double>> w_syn(
        static_cast<size_t>(N) * static_cast<size_t>(K + 1), {0.0, 0.0});
    for (int t = 0; t < N; ++t) {
        for (int q = 0; q < nq; ++q) {
            // The weight of the HARMONIC p = q - 2K (the slot is just an
            // index shift: the harmonic of the slot q is p).
            const double ang = -domega * static_cast<double>(t) *
                               static_cast<double>(q - 2 * K);
            w_dft[static_cast<size_t>(t) * static_cast<size_t>(nq) +
                   static_cast<size_t>(q)] = std::complex<double>(std::cos(ang),
                                                                 std::sin(ang));
        }
        for (int m = 0; m <= K; ++m) {
            const double ang = domega * static_cast<double>(t) *
                               static_cast<double>(m);
            w_syn[static_cast<size_t>(t) * static_cast<size_t>(K + 1) +
                  static_cast<size_t>(m)] = std::complex<double>(std::cos(ang),
                                                                std::sin(ang));
        }
    }

    // ---- the expanded harmonic system ---------------------------------------
    // The unknowns are the complex harmonics −K … K (2K+1 blocks of n) and
    // the system is the block matrix whose (b, m) block is
    //     Ĵ_{p_b − p_m} + j·p_b·ω·Ĉ_{p_b − p_m}
    // with Ĵ_p the DFT of the time-domain dF/dx.
    const int dim = nb * n;
    const int mna_nnz = static_cast<int>(system.jacobianF().values().size());
    const auto& rp = system.jacobianF().rowPointers();
    const auto& ci = system.jacobianF().colIndices();

    // The pattern: the row (b, i) has the columns (m, j) for every MNA
    // entry (i, j) and every block m. The pattern NEVER changes, so the
    // symbolic factorization runs once.
    std::vector<int> row_ptr(static_cast<size_t>(dim) + 1, 0);
    std::vector<int> col_idx;
    std::vector<int> row_off(static_cast<size_t>(dim), 0);
    for (int b = 0; b < nb; ++b) {
        for (int i = 0; i < n; ++i) {
            const int r = b * n + i;
            row_ptr[static_cast<size_t>(r)] = static_cast<int>(col_idx.size());
            for (int p = rp[i]; p < rp[i + 1]; ++p) {
                const int j = ci[static_cast<size_t>(p)];
                for (int m = 0; m < nb; ++m) {
                    col_idx.push_back(m * n + j);
                }
            }
        }
    }
    row_ptr[static_cast<size_t>(dim)] = static_cast<int>(col_idx.size());
    for (int r = 0; r < dim; ++r) {
        row_off[static_cast<size_t>(r)] = row_ptr[static_cast<size_t>(r)];
    }
    const std::vector<int> row_ptr_unperm = row_ptr;

    // ---- the working arrays ------------------------------------------------
    std::vector<double> x_t(static_cast<size_t>(n), 0.0);
    std::vector<double> F_n(static_cast<size_t>(n), 0.0);
    std::vector<double> Q_n(static_cast<size_t>(n), 0.0);
    std::vector<double> B_n(static_cast<size_t>(n), 0.0);
    std::vector<std::complex<double>> R(static_cast<size_t>(dim), {0.0, 0.0});
    std::vector<std::complex<double>> hatG(static_cast<size_t>(nq) *
                                           static_cast<size_t>(mna_nnz));
    std::vector<std::complex<double>> hatC(static_cast<size_t>(nq) *
                                           static_cast<size_t>(mna_nnz));
    std::vector<double> exp_re(col_idx.size(), 0.0);
    std::vector<double> exp_im(col_idx.size(), 0.0);
    std::vector<double> val_re(col_idx.size(), 0.0);
    std::vector<double> val_im(col_idx.size(), 0.0);
    std::vector<double> rhs_re(static_cast<size_t>(dim), 0.0);
    std::vector<double> rhs_im(static_cast<size_t>(dim), 0.0);
    std::vector<double> dx_re(static_cast<size_t>(dim), 0.0);
    std::vector<double> dx_im(static_cast<size_t>(dim), 0.0);
    std::vector<std::complex<double>> X_try(static_cast<size_t>(K + 1) *
                                            static_cast<size_t>(n), {0.0, 0.0});

    // ---- the evaluation ----------------------------------------------------
    // Evaluates the devices at the phase points and transforms the residual
    // back to the frequency domain. `jac` also asks for the DFT of the
    // time-domain Jacobians (the Newton system needs it).
    const auto evaluate = [&](const std::vector<std::complex<double>>& Xc,
                              bool jac) -> double {
        std::fill(R.begin(), R.end(), std::complex<double>(0.0, 0.0));
        if (jac) {
            std::fill(hatG.begin(), hatG.end(), std::complex<double>(0.0, 0.0));
            std::fill(hatC.begin(), hatC.end(), std::complex<double>(0.0, 0.0));
        }

        for (int t = 0; t < N; ++t) {
            // The synthesis: x(t_n) = Re{Σ_m X_m·e^{jmωt_n}}.
            std::fill(x_t.begin(), x_t.end(), 0.0);
            for (int m = 0; m <= K; ++m) {
                const std::complex<double> w =
                    w_syn[static_cast<size_t>(t) * static_cast<size_t>(K + 1) +
                          static_cast<size_t>(m)];
                const double wr = w.real();
                const double wi = w.imag();
                for (int i = 0; i < n; ++i) {
                    const std::complex<double>& xk =
                        Xc[static_cast<size_t>(m) * static_cast<size_t>(n) +
                           static_cast<size_t>(i)];
                    x_t[static_cast<size_t>(i)] += xk.real() * wr - xk.imag() * wi;
                }
            }

            const double time = static_cast<double>(t) * period /
                                static_cast<double>(N);
            const DaeSolution sol(x_t, circuit_, node_map);
            system.evaluateAt(sol, time, F_n, Q_n, B_n);

            // The source excitation at this phase point: the `AC` phasor
            // when the source has one, its waveform otherwise.
            std::fill(B_n.begin(), B_n.end(), 0.0);
            for (const auto& d : circuit_.getDevices()) {
                d->loadHBExcitation(time, omega0, B_n);
            }

            // The DFT of the residual: R_p = (1/N)·Σ_t r(t_t)·e^{-jpωt_t},
            // with r = F + j·p·ω·Q − B.
            for (int q = 0; q < nq; ++q) {
                const int p = q - 2 * K;         // the harmonic of this slot
                if (p < -K || p > K) continue;   // only the solved harmonics
                const std::complex<double> w =
                    w_dft[static_cast<size_t>(t) * static_cast<size_t>(nq) +
                          static_cast<size_t>(q)];
                if (w.real() == 0.0 && w.imag() == 0.0) continue;
                const double jpw = omega0 * p;
                const int b = p + K;              // the block of the harmonic
                for (int i = 0; i < n; ++i) {
                    const double fq = F_n[static_cast<size_t>(i)];
                    const double qq = Q_n[static_cast<size_t>(i)];
                    const double bq = B_n[static_cast<size_t>(i)];
                    R[static_cast<size_t>(b * n + i)] +=
                        w * std::complex<double>(fq - bq, jpw * qq);
                }
            }

            if (!jac) continue;

            // The Jacobians at this phase, and their DFT.
            system.evaluateJacobiansAt(sol, time);
            const auto& gv = system.jacobianF().values();
            const auto& cv = system.jacobianQ().values();
            for (int q = 0; q < nq; ++q) {
                const std::complex<double> w =
                    w_dft[static_cast<size_t>(t) * static_cast<size_t>(nq) +
                          static_cast<size_t>(q)];
                if (w.real() == 0.0 && w.imag() == 0.0) continue;
                std::complex<double>* hg =
                    hatG.data() + static_cast<size_t>(q) * mna_nnz;
                std::complex<double>* hc =
                    hatC.data() + static_cast<size_t>(q) * mna_nnz;
                for (int e = 0; e < mna_nnz; ++e) {
                    hg[e] += w * gv[static_cast<size_t>(e)];
                    hc[e] += w * cv[static_cast<size_t>(e)];
                }
            }
        }

        // The normalization (1/N) of the DFT.
        if (jac) {
            for (auto& v : hatG) v /= static_cast<double>(N);
            for (auto& v : hatC) v /= static_cast<double>(N);
        }

        double res = 0.0;
        for (int b = 0; b < nb; ++b) {
            for (int i = 0; i < n; ++i) {
                res = std::max(res, std::abs(R[static_cast<size_t>(b * n + i)]));
            }
        }
        return res;
    };

    // ---- the Newton iteration ----------------------------------------------
    ComplexSparseLU lu;
    std::vector<int> perm_of_pos;   // permuted position -> original row
    bool ordered = false;
    double res_norm = 0.0;
    int iter = 0;
    bool converged = false;

    for (iter = 0; iter < max_iterations_; ++iter) {
        res_norm = evaluate(X, /*jac=*/true);
        if (std::getenv("EXSPICE_DEBUG_HB")) {
            std::fprintf(stderr, "[HB] iter %d: |R| = %.6g\n", iter + 1, res_norm);
        }
        if (verbose_) {
            std::cout << "[HB] Newton iteration " << (iter + 1)
                      << ": |R| = " << res_norm << std::endl;
        }
        if (res_norm < tolerance_) {
            converged = true;
            break;
        }

        // ---- the expanded Jacobian (in the unpermuted order) --------------
        {
            size_t pos = 0;
            for (int b = 0; b < nb; ++b) {
                const int p_b = b - K;
                for (int i = 0; i < n; ++i) {
                    for (int p = rp[i]; p < rp[i + 1]; ++p) {
                        for (int m = 0; m < nb; ++m) {
                            const int p_m = m - K;
                            const int q = (p_b - p_m) + 2 * K;
                            const std::complex<double>& hg =
                                hatG[static_cast<size_t>(q) * mna_nnz +
                                     static_cast<size_t>(p)];
                            const std::complex<double>& hc =
                                hatC[static_cast<size_t>(q) * mna_nnz +
                                     static_cast<size_t>(p)];
                            const double jpw = omega0 * p_b;
                            const std::complex<double> val =
                                hg + std::complex<double>(0.0, jpw) * hc;
                            exp_re[pos] = val.real();
                            exp_im[pos] = val.imag();
                            ++pos;
                        }
                    }
                }
            }
        }

        // ---- the pivot-safe ordering (decided once, on the real matrix) ----
        if (!ordered) {
            // Priority: the branch rows (the MNA rows after the node rows)
            // make the best partners (their couplings are the ±1 of a branch
            // equation).
            std::vector<int> priority(static_cast<size_t>(dim), 0);
            const int num_node_rows = system.numNodeRows();
            for (int b = 0; b < nb; ++b) {
                for (int i = num_node_rows; i < n; ++i) {
                    priority[static_cast<size_t>(b * n + i)] = 1;
                }
            }
            const StaticRowPermutation perm = StaticRowPermutation::build(
                row_ptr_unperm, col_idx, exp_re, exp_im, priority);
            if (!perm.ok) {
                result.message = "the harmonic system cannot be ordered for "
                                 "factorization: " + perm.reason;
                return result;
            }
            if (verbose_) {
                std::cout << "[HB] pivot-safe ordering: " << perm.weak_diagonals
                          << " column(s) without a self-term reordered ("
                          << perm.swaps << " row swap(s))" << std::endl;
            }
            perm_of_pos = perm.row_of_pos;
            // Apply the permutation to the pattern.
            std::vector<int> new_rp(static_cast<size_t>(dim) + 1, 0);
            std::vector<int> new_ci;
            new_ci.reserve(col_idx.size());
            for (int k = 0; k < dim; ++k) {
                new_rp[static_cast<size_t>(k)] = static_cast<int>(new_ci.size());
                const int src = perm_of_pos[static_cast<size_t>(k)];
                for (int p = row_ptr_unperm[static_cast<size_t>(src)];
                     p < row_ptr_unperm[static_cast<size_t>(src) + 1]; ++p) {
                    new_ci.push_back(col_idx[static_cast<size_t>(p)]);
                }
            }
            new_rp[static_cast<size_t>(dim)] = static_cast<int>(new_ci.size());
            row_ptr.swap(new_rp);
            col_idx.swap(new_ci);
            lu.setPattern(row_ptr, col_idx);
            ordered = true;
        }

        // ---- the values in the permuted order ------------------------------
        {
            size_t dst = 0;
            for (int k = 0; k < dim; ++k) {
                const int r = perm_of_pos[static_cast<size_t>(k)];
                const int off = row_off[static_cast<size_t>(r)];
                const int len = row_ptr_unperm[static_cast<size_t>(r) + 1] - off;
                for (int t = 0; t < len; ++t) {
                    val_re[dst + static_cast<size_t>(t)] = exp_re[static_cast<size_t>(off + t)];
                    val_im[dst + static_cast<size_t>(t)] = exp_im[static_cast<size_t>(off + t)];
                }
                dst += static_cast<size_t>(len);
            }
        }
        if (std::getenv("EXSPICE_DEBUG_HB")) {
            std::fprintf(stderr,
                         "[HB] dim=%d perm[0]=%d perm[1]=%d | A'[0][0]=(%.6g,%.6g) "
                         "A'[0][1]=(%.6g,%.6g) A'[1][0]=(%.6g,%.6g)\n",
                         dim, perm_of_pos[0], perm_of_pos[1],
                         val_re[0], val_im[0],
                         val_re[1], val_im[1],
                         val_re[row_ptr[1]], val_im[row_ptr[1]]);
        }

        // ---- the right-hand side: -R, in the permuted row order ------------
        for (int k = 0; k < dim; ++k) {
            const int r = perm_of_pos[static_cast<size_t>(k)];
            rhs_re[static_cast<size_t>(k)] = -R[static_cast<size_t>(r)].real();
            rhs_im[static_cast<size_t>(k)] = -R[static_cast<size_t>(r)].imag();
        }

        // ---- the solve ------------------------------------------------------
        if (!lu.factorize(0.0, val_re.data(), val_im.data())) {
            result.message = "the harmonic system could not be factorized: " +
                             lu.message();
            return result;
        }
        dx_re = rhs_re;
        dx_im = rhs_im;
        lu.solve(dx_re, dx_im);

        // ---- the step, with a damping when the plain Newton step overshoots -
        // The harmonic balance Newton is notoriously prone to overshooting
        // on a strongly nonlinear circuit, so the step is halved until the
        // residual decreases (at most a few halvings). The acceptance uses
        // the L2 norm of the residual: the max-norm can grow while the
        // solution improves (one badly damped harmonic), and accepting such
        // a step is what makes the iteration diverge.
        //
        // The expanded system's block b carries the harmonic p = b − K, while
        // the stored harmonics are indexed by k = 0 … K (the DC and the
        // positive ones; the negative harmonics are their conjugates and are
        // implied). The block of the harmonic k is b = K + k.
        double alpha = 1.0;
        bool accepted = false;
        for (int attempt = 0; attempt < 10; ++attempt) {
            X_try = X;
            for (int b = 0; b < nb; ++b) {
                const int p = b - K;
                if (p < 0) continue;   // a negative harmonic: implied
                const int k = p;        // the stored harmonic index
                for (int i = 0; i < n; ++i) {
                    const int r = b * n + i;
                    // The step of the block b lives at the permuted position
                    // of the original row r.
                    int kpos = -1;
                    for (int p2 = 0; p2 < dim; ++p2) {
                        if (perm_of_pos[static_cast<size_t>(p2)] == r) {
                            kpos = p2;
                            break;
                        }
                    }
                    if (kpos < 0) continue;
                    X_try[static_cast<size_t>(k) * static_cast<size_t>(n) +
                          static_cast<size_t>(i)] +=
                        alpha * std::complex<double>(
                            dx_re[static_cast<size_t>(kpos)],
                            dx_im[static_cast<size_t>(kpos)]);
                }
            }
            const double res_try = evaluate(X_try, /*jac=*/false);
            if (res_try < res_norm) {
                accepted = true;
                X = X_try;
                if (verbose_ && alpha < 1.0) {
                    std::cout << "[HB]   step damped (alpha = " << alpha
                              << "), |R| = " << res_try << std::endl;
                }
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            // No damping reduced the residual: the iteration is stagnating.
            // The last (tiny) step is kept so the residual norm reported is
            // the one of the accepted point.
            X = X_try;
        }
    }

    result.newton_iterations = iter + 1;
    result.residual = res_norm;
    if (!converged) {
        result.message = "the harmonic balance Newton iteration did not "
                         "converge (|R| = " + std::to_string(res_norm) + ")";
        return result;
    }

    // ---- the results --------------------------------------------------------
    result.harmonics.resize(static_cast<size_t>(K) + 1);
    for (int k = 0; k <= K; ++k) {
        result.harmonics[static_cast<size_t>(k)].resize(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            result.harmonics[static_cast<size_t>(k)][static_cast<size_t>(i)] =
                X[static_cast<size_t>(k) * static_cast<size_t>(n) +
                  static_cast<size_t>(i)];
        }
    }

    // The reconstructed waveform (the real part of the synthesis).
    result.times.resize(static_cast<size_t>(N));
    result.waveform.resize(static_cast<size_t>(N));
    for (int t = 0; t < N; ++t) {
        result.times[static_cast<size_t>(t)] =
            static_cast<double>(t) * period / static_cast<double>(N);
        result.waveform[static_cast<size_t>(t)].assign(static_cast<size_t>(n), 0.0);
        for (int m = 0; m <= K; ++m) {
            const double wr = std::cos(domega * static_cast<double>(t) *
                                       static_cast<double>(m));
            const double wi = std::sin(domega * static_cast<double>(t) *
                                       static_cast<double>(m));
            for (int i = 0; i < n; ++i) {
                const std::complex<double>& xk =
                    X[static_cast<size_t>(m) * static_cast<size_t>(n) +
                      static_cast<size_t>(i)];
                result.waveform[static_cast<size_t>(t)][static_cast<size_t>(i)] +=
                    xk.real() * wr - xk.imag() * wi;
            }
        }
    }

    result.success = true;
    std::ostringstream msg;
    msg << result.newton_iterations << " Newton iteration(s), " << N
        << " time points, " << (K + 1) << " harmonic(s), " << n
        << " unknowns, |R| = " << res_norm;
    result.message = msg.str();
    return result;
}

} // namespace e_XSpice
