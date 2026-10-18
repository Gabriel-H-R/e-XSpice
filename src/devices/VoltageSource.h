/**
 * @author Gabriel Henrique Silva
 * @file VoltageSource.h
 * @brief Independent voltage source — MNA with branch-current unknown.
 *
 * DAE formulation (Xyce-style): a voltage source adds a branch-current
 * unknown i (the current flowing + → − through the source) and a branch
 * equation that enforces  v+ − v− = V(t):
 *
 *     F[pos] += i        F[neg] -= i        F[bra] += (v+ − v−)
 *     B[bra] += V(t)                              (source value, RHS only)
 *     J_F: [pos][bra] = +1   [neg][bra] = -1
 *          [bra][pos] = +1   [bra][neg] = -1
 *
 * The residual of the whole system is F + dQ/dt − B = 0, so the branch row
 * reads  v+ − v− + 0 − V(t) = 0  — exactly the KVL constraint.
 */

#ifndef VOLTAGE_SOURCE_H
#define VOLTAGE_SOURCE_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "SourceWaveform.h"

namespace e_XSpice {

class VoltageSource : public Device {
public:
    VoltageSource(const std::string& name,
                  Node* node_plus,
                  Node* node_minus,
                  SourceWaveform waveform);

    ~VoltageSource() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return 1; }

    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    void seedInitialGuess(std::vector<double>& x) const override;

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    void loadF(const DaeSolution& x, double time,
               std::vector<double>& F) const override;

    void loadB(double time, std::vector<double>& B) const override;

    /**
     * @brief Small-signal phasor on the branch row (AC analysis).
     *        The KVL constraint itself is already in dF/dx, so only the
     *        source phasor (mag·e^{j·phase}) enters the right-hand side.
     */
    void loadACExcitation(double omega, std::vector<double>& Breal,
                          std::vector<double>& Bimag) const override;

    void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    void requestBreakpoints(std::vector<double>& times, double tstop) const override;

    double maxTimeStepHint(double t) const override;

    /**
     * @brief Recomputes cached outputs from the global solution (exact on
     *        every rank, including ranks that do not own this device).
     */
    void refreshOutputs(const DaeSolution& x, double time) const override;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    bool   isLinear() const override { return true; }
    double getCurrent() const override { return last_current_; }  // + → − through source
    double getPower()   const override;                            // delivered (V·i)

    // Getters / setters
    double getVoltage() const { return waveform_.value(last_time_); }
    const SourceWaveform& waveform() const { return waveform_; }
    void setWaveform(const SourceWaveform& w) { waveform_ = w; }

private:
    SourceWaveform waveform_;
    mutable double last_current_ = 0.0;  // last accepted branch current (A)
    mutable double last_vcap_   = 0.0;
    mutable double last_time_   = 0.0;
};

} // namespace e_XSpice

#endif // VOLTAGE_SOURCE_H
