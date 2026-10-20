/**
 * @author Gabriel Henrique Silva
 * @file TransientSensitivityAnalysis.h
 * @brief Transient sensitivity analysis (`.SENS TRAN`): dO/dp(t) with
 *        respect to circuit parameters p, by the adjoint method.
 *
 * Card
 * ----
 *   .SENS TRAN V(out) R1,R2
 *   .SENS TRAN objfunc=V(out) param=R1,R2
 *
 * Method (adjoint, step by step)
 * -----------------------------
 * At every accepted time step the transient engine solves the nonlinear DAE
 *
 *      F(x,p) + s·(Q(x,p) − Q_n) = B(p)
 *
 * with the effective Jacobian `J = dF/dx + s·dQ/dx`. Differentiating the
 * CONVERGED residual with respect to p gives
 *
 *      J·(dx/dp) = w ,   w = ∂B/∂p − (∂G/∂p)·X − s·(∂C/∂p)·X
 *
 * and, for an objective `O(X)`,
 *
 *      dO/dp = λ^T·w ,   λ = J^{-T}·(∂O/∂X)
 *
 * so ONE adjoint solve per time step serves EVERY parameter: the cost of the
 * whole sweep is independent of the number of parameters (the direct method
 * would need one solve per parameter per step). The perturbation `w` comes
 * from the device layer, which is what makes the engine device-agnostic and
 * what makes the DC, AC and transient sensitivity analyses the SAME
 * computation with a different integration coefficient.
 *
 * The analysis hooks into the transient engine through its step observer
 * instead of re-implementing the time loop, so the reported sensitivities
 * belong to EXACTLY the same trajectory a plain `.TRAN` run produces
 * (breakpoints, step rejection and the LTE controller included).
 *
 * Method (direct, with the recursion of the trajectory)
 * ------------------------------------------------------
 * The trajectory is a recursion `x_{n+1} = Φ(x_n, p)`, whose sensitivity
 * obeys `Y_{n+1} = ∂Φ/∂x·Y_n + ∂Φ/∂p`. The objective at step n is
 * `O(x_n)`, so `dO/dp = (∂O/∂x)·Y_n` and the sensitivity VECTOR has to be
 * propagated: one solve per parameter per step (the adjoint of a recursion
 * is a backward sweep, which cannot produce the sensitivity of every time
 * point of a forward simulation — the adjoint is the right tool for a single
 * objective at the END of the trajectory).
 *
 *      J_n·Y_n = w_n
 *
 * with `J_n = dF/dx + s·dQ/dx` (the effective Jacobian of the step) and the
 * right-hand-side perturbation
 *
 *      w_n = ∂B/∂p − (∂G/∂p)·X_n − s·(∂C/∂p)·X_n   (the device's own part)
 *            + s·(dQ/dx)·Y_{n-1}                     (the charge history)
 *            + e_{n-1}                              (trapezoidal rule)
 *
 * The perturbation of the device layer
 * (`Device::loadSensitivity`) is the SAME mathematical object as in the DC
 * and AC flavours, which is what makes the three analyses one computation
 * with a different integration coefficient.
 *
 * Covered parameters (analytical w): R, C, L, V, I sources.
 * Nonlinear devices are not supported yet.
 */

#ifndef TRANSIENT_SENSITIVITY_ANALYSIS_H
#define TRANSIENT_SENSITIVITY_ANALYSIS_H

#include "../core/Circuit.h"
#include "TransientAnalysis.h"

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace e_XSpice {

class DaeSystem;

struct TransientSensitivityResult {
    bool success = false;
    std::string message;

    std::string output_var;
    std::vector<std::string> param_names;   // "R1:R", ...

    struct Point {
        double time = 0.0;
        double obj_value = 0.0;              // O(t)
        std::vector<double> sens;            // dO/dp(t), aligned with param_names
        bool valid = true;
    };
    std::vector<Point> points;
};

class TransientSensitivityAnalysis {
public:
    explicit TransientSensitivityAnalysis(Circuit& circuit);
    ~TransientSensitivityAnalysis();

    void setVariation(double tstart, double tstop, double tstep) {
        tstart_ = tstart; tstop_ = tstop; tstep_ = tstep;
    }
    void setOutputVar(const std::string& var) { output_var_ = var; }
    void addParameter(const std::string& device, const std::string& param) {
        params_.emplace_back(device, param);
    }
    void setMaxIterations(int n) { max_iterations_ = n; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool v) { verbose_ = v; }
    void setCsvStream(std::ostream* os) { csv_ = os; }
    void setIntegrationMethod(IntegrationMethod m) { method_ = m; }
    TransientSensitivityResult run(bool use_mpi = true);

private:
    Circuit& circuit_;
    double tstart_ = 0.0;
    double tstop_ = 1e-3;
    double tstep_ = 1e-4;

    std::string output_var_;
    std::vector<std::pair<std::string, std::string>> params_;

    int    max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool   verbose_ = false;
    IntegrationMethod method_ = IntegrationMethod::TRAPEZOIDAL;

    std::ostream* csv_ = nullptr;
    bool csv_header_written_ = false;
};

} // namespace e_XSpice

#endif // TRANSIENT_SENSITIVITY_ANALYSIS_H
