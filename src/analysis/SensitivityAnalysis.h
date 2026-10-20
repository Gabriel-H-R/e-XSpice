/**
 * @author Gabriel Henrique Silva
 * @file SensitivityAnalysis.h
 * @brief DC sensitivity analysis (`.SENS`): dO/dp at the operating point,
 *        by the adjoint method.
 *
 * Card
 * ----
 *   .SENS <outvar> [param1, param2, ...]          (ngspice style)
 *   .SENS objfunc=<outvar> param=<p1,p2,...>      (Xyce style)
 *
 *   - `outvar` is `V(node)`, `V(n1,n2)` or `I(device)`.
 *   - params are device names (`R1`) or `device:param` (`R1:R`). If omitted,
 *     every device's sensitivity parameter is used.
 *
 * Method (adjoint)
 * ----------------
 * At the operating point the small-signal system is `A·X = B` with `A = dF/dx`
 * (the static Jacobian). For an objective `O(X)`:
 *
 *     λ = A^{-T}·(∂O/∂X)                       (ONE solve per objective)
 *     dO/dp = λ^T·w ,   w = ∂B/∂p − (∂A/∂p)·X  (per parameter, no solve)
 *
 * so the cost is `O(1 solve + n_params·eval)` — INDEPENDENT of the number of
 * parameters. Xyce supports both the direct method (`O(n_params·solve)`) and
 * the adjoint one; this engine always uses the adjoint, which is the scalable
 * choice and the reason the whole sweep of a large parameter set costs the
 * same as a single parameter.
 *
 * The perturbation `w` is produced by the device layer
 * (`Device::loadSensitivity`), so the analysis engine is device-agnostic: a
 * new device becomes differentiable just by implementing that one hook.
 *
 * Covered parameters today (analytical w):
 *   - Resistor: `R`
 *   - Capacitor: `C` (no effect at DC — a capacitor is an open circuit)
 *   - Inductor: `L` (no effect at DC — an inductor is a short circuit)
 *   - Voltage source: `DC`/`V`
 *   - Current source: `DC`/`I`
 * Nonlinear devices (diodes, BJTs, MOSFETs) are not supported yet: their
 * residual is not linear in the unknowns, so the device-level derivative is
 * not the matrix-level one.
 */

#ifndef SENSITIVITY_ANALYSIS_H
#define SENSITIVITY_ANALYSIS_H

#include "../core/Circuit.h"

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace e_XSpice {

struct SensitivityResult {
    bool success = false;
    std::string message;

    std::string output_var;
    std::vector<std::string> param_names;   // "R1:R", ...
    std::vector<double> sensitivities;     // dO/dp, aligned with param_names
};

class SensitivityAnalysis {
public:
    explicit SensitivityAnalysis(Circuit& circuit);
    ~SensitivityAnalysis();

    void setOutputVar(const std::string& var) { output_var_ = var; }
    void addParameter(const std::string& device, const std::string& param) {
        params_.emplace_back(device, param);
    }
    void setVerbose(bool v) { verbose_ = v; }
    void setCsvStream(std::ostream* os) { csv_ = os; }

    SensitivityResult run();

private:
    Circuit& circuit_;
    std::string output_var_;
    std::vector<std::pair<std::string, std::string>> params_;
    bool verbose_ = false;
    std::ostream* csv_ = nullptr;
};

} // namespace e_XSpice

#endif // SENSITIVITY_ANALYSIS_H
