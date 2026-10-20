/**
 * @author Gabriel Henrique Silva
 * @file SensitivityContext.h
 * @brief Context of a sensitivity evaluation (`.SENS`).
 *
 * The small-signal sensitivity of an output O with respect to a circuit
 * parameter p is ALWAYS the same mathematical object, whatever the analysis
 * that produced the operating point:
 *
 *      A(ξ) · X = B(ξ)                (ξ = 0 for DC, j·w for AC, s for
 *                                       the transient integrator)
 *      dO/dp = λ^H · w                (λ = A^{-H}·(∂O/∂X))
 *
 * with the RIGHT-HAND-SIDE PERTURBATION
 *
 *      w = ∂B/∂p − (∂A/∂p)·X
 *
 * and, for the transient (where the effective Jacobian is J = G + s·C and the
 * residual is R = F + s·Q − B):
 *
 *      w = ∂B/∂p − (∂G/∂p)·X − s·(∂C/∂p)·X
 *
 * The device layer produces `w` (one vector, in the DAE row space) and the
 * analysis engines only have to (a) build the right A, (b) solve the adjoint
 * system once and (c) take the scalar product. That is why the SAME device
 * hook serves the DC, AC and transient sensitivity analyses: only the value
 * of the integration coefficient ξ changes.
 *
 * The hook lives on Device (`loadSensitivity`) and this struct tells the
 * device which flavour of small-signal problem it is answering:
 *
 *   - DC        : ξ = 0, the perturbation is the static one.
 *   - AC        : ξ = j·w, the solution X is a PHASOR (complex).
 *   - TRANSIENT : ξ = 2/dt (trapezoidal) or 1/dt (backward Euler), X is real
 *                 and the analysis re-evaluates it at every accepted step.
 */

#ifndef SENSITIVITY_CONTEXT_H
#define SENSITIVITY_CONTEXT_H

#include <complex>

namespace e_XSpice {

/**
 * @brief What the sensitivity evaluation is being asked for.
 */
struct SensitivityContext {
    enum class Kind {
        DC,          ///< Operating point: ξ = 0.
        AC,          ///< Small signal: ξ = j·ω, X is a phasor.
        TRANSIENT    ///< Time domain: ξ = 2/dt or 1/dt, X is real.
    };

    Kind kind = Kind::DC;

    /** Simulation time (s). DC: 0. AC: 0. TRANSIENT: the step time. */
    double time = 0.0;

    /** Angular frequency 2·π·f (AC only). */
    double omega = 0.0;

    /**
     * @brief Integration coefficient (transient only): 2/dt for the
     *        trapezoidal rule, 1/dt for backward Euler.
     */
    double s = 0.0;

    /**
     * @brief The integration coefficient as a complex number.
     *
     * DC → 0, AC → j·ω, TRANSIENT → s (real).
     */
    std::complex<double> coeff() const {
        switch (kind) {
            case Kind::AC:        return std::complex<double>(0.0, omega);
            case Kind::TRANSIENT: return std::complex<double>(s, 0.0);
            case Kind::DC:        break;
        }
        return std::complex<double>(0.0, 0.0);
    }

    /** True when the solution the device must read is a phasor. */
    bool isPhasor() const { return kind == Kind::AC; }
};

} // namespace e_XSpice

#endif // SENSITIVITY_CONTEXT_H
