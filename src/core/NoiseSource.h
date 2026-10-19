/**
 * @author Gabriel Henrique Silva
 * @file NoiseSource.h
 * @brief Noise-generator contract: what a device reports about its own noise.
 *
 * ================================================================
 * WHY THIS IS A REPORT, NOT A STAMP
 * ================================================================
 * Noise does not enter the linear system. The frequency-domain small-signal
 * system is solved ONCE per frequency, and every noise generator of the
 * circuit is a CURRENT SOURCE in parallel with its device whose power
 * spectral density is known at the operating point. The analysis then needs,
 * for each of them, the transfer from its two terminals to the output:
 *
 *      V_out(f) = sum_i  H_i(f) · I_i(f)
 *
 * with H_i = A^-1·(e_pos - e_neg). Instead of solving one linear system per
 * generator (N solves per frequency), the ADJOINT is solved once:
 *
 *      A^T · x* = e_out          =>     H_i = (x*[pos] - x*[neg])
 *
 * by symmetry of A, at the cost of a single triangular solve. That is the
 * whole reason a noise sweep scales with the circuit size and not with
 * (sources x unknowns): the per-generator cost is a dot product.
 *
 * ================================================================
 * THE MODEL: S(f) = psd_a2 · f^(-f_exponent)
 * ================================================================
 * Every standard generator fits this two-parameter form:
 *   - thermal (Nyquist)      S = 4kT/R              -> exponent 0
 *   - shot                    S = 2q·I               -> exponent 0
 *   - flicker (1/f, 1/f^NF)   S = KF·I^AF/(Cox·Leff^2·f^NF) -> exponent NF
 * so the frequency-dependent part is computed ONCE PER FREQUENCY as a single
 * multiply, with no per-device virtual call inside the sweep (Xyce calls
 * getNoiseSources() on every device at every frequency).
 *
 * Sign/geometry convention: the generator is a Norton pair injecting +i into
 * `pos` and -i into `neg`, i.e. the current flows from `pos` to `neg` THROUGH
 * the device, exactly like the resistor stamp. `pos`/`neg` are DAE row indices
 * (Device::rowOf()), and the analysis maps them onto the ground-eliminated
 * small-signal system; the ground node is mapped to -1 by the analysis.
 */

#ifndef NOISE_SOURCE_H
#define NOISE_SOURCE_H

#include <string>
#include <vector>

namespace e_XSpice {

class DaeSolution;

/**
 * @brief One noise generator of one device, evaluated at the operating point.
 */
struct NoiseSource {
    /** Generator label, e.g. "thermal", "shot", "flicker", "rb", "ic".
     *  Becomes part of the CSV column name together with the device name. */
    std::string label;

    /** DAE row of the + terminal (rowOf(node[0])); -1 for ground. */
    int pos = -1;

    /** DAE row of the - terminal (rowOf(node[1])); -1 for ground. */
    int neg = -1;

    /** Current-noise prefactor S(1 Hz) in A^2/Hz (SPICE's "noise density"). */
    double psd_a2 = 0.0;

    /** S(f) = psd_a2 * f^(-f_exponent). 0 for white noise (thermal, shot). */
    double f_exponent = 0.0;
};

/**
 * @brief Read-only context handed to Device::loadNoiseSources().
 */
struct NoiseContext {
    /** Circuit temperature (K) — the noise is thermal after all. */
    double temperature = 300.15;
    /** Boltzmann constant (J/K). */
    double boltzmann = 1.380649e-23;
    /** Elementary charge (C). */
    double charge = 1.602176634e-19;
    /** Converged operating point (node voltages + branch currents). */
    const DaeSolution* op = nullptr;
};

} // namespace e_XSpice

#endif // NOISE_SOURCE_H