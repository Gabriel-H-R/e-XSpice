/**
 * @author Gabriel Henrique Silva
 * @file HBAnalysis.h
 * @brief Harmonic Balance analysis (`.HB`): the steady-state periodic
 *        response of a (nonlinear) circuit, in the frequency domain.
 *
 * Card
 * ----
 *   .HB <fundamental frequency>
 *   .OPTIONS HBINT NUMFREQ=<K>      (number of harmonics, default 3)
 *
 * Method
 * ------
 * The circuit is described by the same DAE as every other analysis,
 *
 *      C(x)·ẋ + F(x) = b(t)
 *
 * and the periodic steady state is written as a Fourier series
 *
 *      x(t) = Σ_k X_k·e^{jkωt}        (k = 0 … K)
 *
 * The unknowns are therefore the HARMONIC COEFFICIENTS X_k of the MNA
 * solution, and the equations are the Fourier coefficients of the residual:
 *
 *      R_k = F_k + jkω·Q_k − B_k = 0        (k = 0 … K)
 *
 * where F_k, Q_k and B_k are the discrete Fourier transforms of the
 * device residual, charge and source vectors evaluated at the phase points
 * of one period:
 *
 *      x(t_n) = Σ_m X_m·e^{jmωt_n}         (the synthesis, n = 0 … N-1)
 *      F_k    = (1/N)·Σ_n F(x(t_n))·e^{−jkωt_n}      (the analysis)
 *
 * The Jacobian of that system is the CONVOLUTION of the time-domain
 * Jacobian with the DFT:
 *
 *      ∂R_k/∂X_m = (1/N)·Σ_n [∂F/∂x + jkω·∂Q/∂x](t_n)·e^{j(m−k)ωt_n}
 *
 * which is what the engine builds: the per-entry DFT of the time-domain
 * Jacobian is computed ONCE per Newton iteration and the blocks are then
 * read off it — the cost is O(nnz·N·K) for the DFTs plus O(K²·nnz) for the
 * assembly, instead of the O(N·K²·nnz) a direct accumulation would cost.
 *
 * Why this engine is not a copy of Xyce's
 * ---------------------------------------
 *  - It is DEVICE-AGNOSTIC: the HB reuses the very same device contract as
 *    the transient and the AC analyses (F, Q, B and the two Jacobians), so
 *    every device that works in the time domain works here — no per-device
 *    HB path, no device-specific restriction list.
 *  - The Newton limiting (PNJLIM) is a convergence aid of the time-domain
 *    solver and it CHANGES the residual, so the HB turns it off (Xyce has
 *    the same switch, `.OPTIONS HBINT VOLTLIM`).
 *  - The initial guess does not need a transient run: the operating point
 *    gives the DC harmonic and the AC solution at the fundamental gives
 *    the first harmonic (Xyce's default is a full transient-assisted run,
 *    which costs a whole transient simulation just to seed the Newton).
 *  - The harmonic system is assembled on the DAE's FIXED sparsity pattern,
 *    so the symbolic factorization runs once and every Newton iteration only
 *    refills the values.
 *
 * Limitations (documented, not silent):
 *  - Single tone. The multi-tone case (several fundamentals, intermodulation
 *    products) needs a non-commensurate sampling grid and is not implemented.
 *  - The sources must be periodic: a DC value (the DC harmonic) plus the
 *    `AC` phasor (the fundamental). A time-dependent waveform (PULSE, SIN)
 *    is not periodic and is not supported.
 *  - The Newton solve is a dense block solve (the expanded harmonic system),
 *    so the reachable problem size is the one a direct solver can take; the
 *    Krylov alternative for large systems is the natural next step.
 */

#ifndef HB_ANALYSIS_H
#define HB_ANALYSIS_H

#include "../core/Circuit.h"

#include <complex>
#include <iosfwd>
#include <string>
#include <vector>

namespace e_XSpice {

struct HBResult {
    bool success = false;
    std::string message;

    double fundamental = 0.0;      ///< fundamental frequency (Hz)
    int    num_harmonics = 0;      ///< K (the harmonics 0 … K are computed)
    int    num_time_points = 0;    ///< N (the samples per period)
    int    newton_iterations = 0;
    double residual = 0.0;         ///< final residual norm

    /**
     * @brief The harmonic coefficients of the MNA solution.
     *
     * harmonics[k][i] is the coefficient X_k of the i-th unknown, for
     * k = 0 … K (the DC harmonic is k = 0).
     */
    std::vector<std::vector<std::complex<double>>> harmonics;

    /// The phase points of the period (s) and the reconstructed waveform.
    std::vector<double> times;
    std::vector<std::vector<double>> waveform;   // [n][unknown]

    /// The DAE layout, so the caller can map the unknowns to node names.
    std::vector<int> node_map;        ///< node id -> DAE row
    int num_unknowns = 0;
};

class HBAnalysis {
public:
    explicit HBAnalysis(Circuit& circuit);
    ~HBAnalysis();

    void setFundamental(double f) { fundamental_ = f; }
    void setHarmonics(int k) { harmonics_ = k; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setMaxIterations(int n) { max_iterations_ = n; }
    void setMaxTimePoints(int n) { max_time_points_ = n; }
    void setVerbose(bool v) { verbose_ = v; }
    void setCsvStream(std::ostream* os) { csv_ = os; }

    HBResult run(bool use_mpi = true);

private:
    Circuit& circuit_;
    double fundamental_ = 1e3;   ///< Hz
    int    harmonics_ = 3;       ///< K
    double tolerance_ = 1e-6;
    int    max_iterations_ = 100;
    int    max_time_points_ = 0; ///< 0 = automatic
    bool   verbose_ = false;
    std::ostream* csv_ = nullptr;
};

} // namespace e_XSpice

#endif // HB_ANALYSIS_H
