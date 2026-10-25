/**
 * @author Gabriel Henique Silva
 * @file CurrentSource.cpp
 * @brief Independent current source — pure RHS injection (B vector).
 */

#include "CurrentSource.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"

namespace e_XSpice {

CurrentSource::CurrentSource(const std::string& name,
                             Node* node_plus,
                             Node* node_minus,
                             SourceWaveform waveform)
    : Device(name, DeviceType::CURRENT_SOURCE)
    , waveform_(std::move(waveform))
{
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    if (node_plus)  node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void CurrentSource::registerJacobianPattern(std::vector<int>& rows,
                                            std::vector<int>& cols) const {
    (void)rows; (void)cols;   // no matrix entries at all
}

void CurrentSource::loadB(double time, std::vector<double>& B) const {
    last_time_ = time;
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const double i = waveform_.value(time);

    if (ip >= 0) B[ip] -= i;
    if (in >= 0) B[in] += i;
}

void CurrentSource::loadACExcitation(double omega,
                                    std::vector<double>& Breal,
                                    std::vector<double>& Bimag) const {
    (void)omega;  // a source phasor is frequency independent by definition
    if (!waveform_.hasAC()) return;       // silent source: no AC excitation
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    double re = 0.0, im = 0.0;
    waveform_.acPhasor(re, im);
    // Same sign convention as loadB: the current leaves the + node inside
    // the source, so it enters the KCL rows as -I at pos and +I at neg.
    if (ip >= 0) {
        Breal[static_cast<size_t>(ip)] -= re;
        Bimag[static_cast<size_t>(ip)] -= im;
    }
    if (in >= 0) {
        Breal[static_cast<size_t>(in)] += re;
        Bimag[static_cast<size_t>(in)] += im;
    }
}

void CurrentSource::requestBreakpoints(std::vector<double>& times,
                                       double tstop) const {
    waveform_.breakpoints(times, tstop);
}

double CurrentSource::maxTimeStepHint(double t) const {
    return waveform_.maxTimeStep(t);
}

double CurrentSource::getPower() const {
    const double vp = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double vn = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    // The current flows from node_plus to node_minus INSIDE the source, so
    // the power delivered to the circuit is −I·(v+ − v−).
    return -waveform_.value(last_time_) * (vp - vn);
}

// ---------------------------------------------------------------
// Sensitivity (`.SENS`): d(F-B)/dI = +1 at pos, -1 at neg.
// ---------------------------------------------------------------
// Sensitivity (`.SENS`): the right-hand-side perturbation w.
//
// The current source has NO static residual (F is empty) and enters the
// system only through B:
//
//   B[pos] -= I,  B[neg] += I
//
// so w = dB/dI is -1 at pos and +1 at neg (DC / transient). In the AC
// flavour the parameter is the AC magnitude: B_ac carries -I_ac at pos,
// hence dB_ac/d(mag) = -e^{jφ} at pos.
// ---------------------------------------------------------------
// ---------------------------------------------------------------
// Harmonic balance excitation.
//
// A source carries two ALTERNATIVE excitations: the time-domain waveform
// (what the transient uses) and the `AC` phasor (what the frequency-domain
// analyses use). The HB is a frequency-domain analysis, so a source with an
// `AC` spec is driven by the phasor: its time-domain contribution is
// Re{B·e^{jωt}}, and the DFT of that gives exactly the phasor back.
// ---------------------------------------------------------------
void CurrentSource::loadHBExcitation(double time, double omega,
                                     std::vector<double>& B) const {
    if (!waveform_.hasAC()) {
        loadB(time, B);
        return;
    }
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    double re = 0.0, im = 0.0;
    waveform_.acPhasor(re, im);
    const double c = std::cos(omega * time);
    const double s = std::sin(omega * time);
    const double val = re * c - im * s;
    if (ip >= 0) B[static_cast<size_t>(ip)] -= val;
    if (in >= 0) B[static_cast<size_t>(in)] += val;
}

bool CurrentSource::loadSensitivity(const std::string& param_name,
                                    const DaeSolution& x,
                                    const SensitivityContext& ctx,
                                    std::vector<double>& w_re,
                                    std::vector<double>& w_im,
                                    std::vector<double>& dqdp_re) const {
    (void)x; (void)dqdp_re;   // a source has no state
    if (param_name != "DC" && param_name != "I" && param_name != "dc" &&
        param_name != "i") {
        return false;
    }
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);

    double dre = -1.0;
    double dim = 0.0;
    if (ctx.kind == SensitivityContext::Kind::AC) {
        if (!waveform_.hasAC()) return false;   // silent source: no AC gain
        const double phi = waveform_.acPhaseRad();
        dre = -std::cos(phi);
        dim = -std::sin(phi);
    }
    if (ip >= 0) {
        w_re[static_cast<size_t>(ip)] += dre;
        w_im[static_cast<size_t>(ip)] += dim;
    }
    if (in >= 0) {
        w_re[static_cast<size_t>(in)] -= dre;
        w_im[static_cast<size_t>(in)] -= dim;
    }
    return true;
}

} // namespace e_XSpice
