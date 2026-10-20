/**
 * @author Gabriel Henrique Silva
 * @file CurrentSource.h
 * @brief Independent current source — pure RHS injection (B vector).
 *
 * An independent current source adds NO unknowns and NO matrix entries:
 * the entire contribution is the injection  I(t)  in the residual RHS
 * (Xyce: B[pos] −= I ; B[neg] += I).
 *
 * Convention: positive current flows from node_plus to node_minus
 * THROUGH the source, i.e. it is injected into node_minus and pulled
 * out of node_plus (SPICE convention).
 */

#ifndef CURRENT_SOURCE_H
#define CURRENT_SOURCE_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "SourceWaveform.h"

namespace e_XSpice {

class CurrentSource : public Device {
public:
    CurrentSource(const std::string& name,
                  Node* node_plus,
                  Node* node_minus,
                  SourceWaveform waveform);

    ~CurrentSource() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return 0; }

    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    void requestBreakpoints(std::vector<double>& times, double tstop) const override;

    double maxTimeStepHint(double t) const override;

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    void loadB(double time, std::vector<double>& B) const override;

    /**
     * @brief Small-signal phasor injected into the two nodes (AC analysis).
     *        A current source has no Jacobian contribution at all: its whole
     *        effect on the small-signal solution is the right-hand side.
     */
    void loadACExcitation(double omega, std::vector<double>& Breal,
                          std::vector<double>& Bimag) const override;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    bool   isLinear() const override { return true; }
    double getCurrent() const override { return waveform_.value(last_time_); }
    double getPower()   const override;   // V·I delivered by the source

    // Sensitivity (`.SENS`)
    std::vector<std::string> sensitivityParameters() const override {
        return {"DC", "I"};
    }
    bool loadSensitivity(const std::string& param_name, const DaeSolution& x,
                         const SensitivityContext& ctx,
                         std::vector<double>& w_re,
                         std::vector<double>& w_im,
                         std::vector<double>& dqdp_re) const override;

    // Getters / setters
    double getValue(double t = 0.0) const { return waveform_.value(t); }
    const SourceWaveform& waveform() const { return waveform_; }
    void setWaveform(const SourceWaveform& w) { waveform_ = w; }

private:
    SourceWaveform waveform_;
    mutable double last_time_ = 0.0;
};

} // namespace e_XSpice

#endif // CURRENT_SOURCE_H
