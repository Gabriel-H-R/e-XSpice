/**
 * @author Gabriel Henrique Silva
 * @file ACSensitivityAnalysis.h
 * @brief AC small-signal sensitivity (`.SENS` with `.AC`): dO/dp(ω) for an
 *        output expression O w.r.t. circuit parameters p, by the adjoint
 *        method, over a frequency sweep.
 *
 * Card
 * ----
 *   .AC DEC 10 1 1e6
 *   .SENS AC V(out) R1,R2
 *
 *   (or Xyce style: .SENS acobjfunc=V(out) param=R1,R2)
 *
 * Method
 * ------
 * At each frequency ω the small-signal system is A(ω) = G + jωC. For an
 * objective O(X(ω)):
 *
 *     λ(ω) = A(ω)^{-H}·(∂O/∂X)          (ONE adjoint solve per objective)
 *     dO/dp(ω) = λ(ω)^H·w ,  w = ∂B/∂p − (∂A/∂p)·X
 *
 * Cost: O(n_freqs × (1 factorization + n_objfuncs × solve)), INDEPENDENT of
 * the number of parameters — the reason the adjoint is used instead of the
 * direct method (one solve per parameter per frequency).
 *
 * The perturbation `w` comes from the device layer
 * (`Device::loadSensitivity`), so the engine is device-agnostic and the DC,
 * AC and transient sensitivity analyses share the very same formulation.
 *
 * Covered parameters (analytical w): R, C, L, V, I sources.
 * Nonlinear devices are not supported yet.
 */

#ifndef AC_SENSITIVITY_ANALYSIS_H
#define AC_SENSITIVITY_ANALYSIS_H

#include "../core/Circuit.h"
#include "ACAnalysis.h"

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace e_XSpice {

struct ACSensitivityResult {
    bool success = false;
    std::string message;
    std::string variation = "DEC";
    std::string output_var;
    std::vector<std::string> param_names;
    // Per frequency: objective magnitude/phase, then per-param sensitivity
    // magnitude/phase.
    struct Point {
        double frequency = 0.0;
        double obj_mag = 0.0;
        double obj_phase_deg = 0.0;
        std::vector<double> sens_mag;
        std::vector<double> sens_phase_deg;
        bool valid = true;
    };
    std::vector<Point> points;
};

class ACSensitivityAnalysis {
public:
    explicit ACSensitivityAnalysis(Circuit& circuit);
    ~ACSensitivityAnalysis();

    void setVariation(ACVariation v) { variation_ = v; }
    void setPoints(int n) { points_ = n; }
    void setFrequencyRange(double fstart, double fstop) {
        fstart_ = fstart; fstop_ = fstop;
    }
    void setOutputVar(const std::string& var) { output_var_ = var; }
    void addParameter(const std::string& device, const std::string& param) {
        params_.emplace_back(device, param);
    }
    void setMaxIterations(int n) { max_iterations_ = n; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool v) { verbose_ = v; }
    void setPointDistribution(bool on) { point_distribution_ = on; }
    void setMaxRetainedPoints(size_t n) { max_retained_ = n; }
    void setCsvStream(std::ostream* os) { csv_ = os; }

    ACSensitivityResult run(bool use_mpi = true);

private:
    Circuit& circuit_;
    ACVariation variation_ = ACVariation::DEC;
    int    points_ = 10;
    double fstart_ = 1.0;
    double fstop_  = 1e6;

    std::string output_var_;
    std::vector<std::pair<std::string, std::string>> params_;

    int    max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool   verbose_ = false;
    bool   point_distribution_ = false;
    size_t max_retained_ = 0;

    std::ostream* csv_ = nullptr;
    bool csv_header_written_ = false;

    static constexpr int kMaxFrequencies_ = 1000000;
};

} // namespace e_XSpice

#endif // AC_SENSITIVITY_ANALYSIS_H
