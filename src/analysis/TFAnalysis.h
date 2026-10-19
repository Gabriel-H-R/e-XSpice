/**
 * @author Gabriel Henrique Silva
 * @file TFAnalysis.h
 * @brief DC transfer-function analysis (`.TF`): the DC small-signal gain, the
 *        input resistance and the output resistance, computed from a single
 *        linearization of the operating point.
 *
 * Card
 * ----
 *     .TF <outvar> <insrc>
 *
 *   - `outvar` is the small-signal output: `V(node)`, `V(n1,n2)` or
 *     `I(device)` (a voltage source, an inductor, or a resistor).
 *   - `insrc` names the independent source the input is applied to.
 *
 * Method
 * ------
 * At the DC operating point the small-signal system is `A = dF/dx` (the
 * Jacobian, with the reactive part set to zero — i.e. `ω = 0`). Three solves
 * on the *same* factorization give everything:
 *
 *   1. drive the input source with a unit phasor and solve `A·x = b_in`:
 *        gain = x[out]  (the input source is 1 by construction, so the ratio
 *        is just the output variable)
 *        Rin  = V(src) / I(src)  — for a voltage source `1 / x[branch]`,
 *        for a current source `x[pos] − x[neg]`.
 *   2. inject a unit current at the output node (with the input source in its
 *      small-signal state — a short for a voltage source, an open for a
 *      current source, which is already what `A` encodes) and solve
 *      `A·x = b_out`:
 *        Rout = x[out_pos] − x[out_neg].
 *
 * The input source is driven with a unit phasor regardless of its `AC` spec,
 * so `.TF` works on circuits that carry no AC information at all — the same
 * convention as the noise analysis.
 *
 * This reuses `SmallSignalLinearization` (operating point + ground elimination
 * + pivot-safe ordering) and the complex LU factorized once at `ω = 0`, so the
 * cost is `O(N)` for the OP plus a handful of triangular solves. No sweep, no
 * frequency loop, no parallelism needed.
 *
 * Note on Xyce: Xyce has no dedicated `.TF` analysis; this is implemented as
 * the standard SPICE3/NGSPICE DC small-signal transfer function.
 */

#ifndef TF_ANALYSIS_H
#define TF_ANALYSIS_H

#include "../core/Circuit.h"

#include <string>

namespace e_XSpice {

/**
 * @brief Result of a `.TF` analysis: a handful of scalars.
 */
struct TFResult {
    bool success = false;
    std::string message;

    std::string output_var;     // as given, e.g. "V(out)" or "I(V1)"
    std::string input_source;  // as given, e.g. "VIN"

    double transfer_function = 0.0;   // output / input (dimensionless or V/A or A/V)
    double input_resistance  = 0.0;   // ohms
    double output_resistance = 0.0;   // ohms
};

class TFAnalysis {
public:
    explicit TFAnalysis(Circuit& circuit);
    ~TFAnalysis();

    void setOutputVar(const std::string& var) { output_var_ = var; }
    void setInputSource(const std::string& src) { input_source_ = src; }

    void setMaxIterations(int n) { max_iterations_ = n; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool v) { verbose_ = v; }

    TFResult run();

private:
    Circuit& circuit_;
    std::string output_var_;
    std::string input_source_;

    int    max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool   verbose_ = false;
};

} // namespace e_XSpice

#endif // TF_ANALYSIS_H