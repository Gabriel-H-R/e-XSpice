/**
 * @author Gabriel Henique Silva
 * @file SourceWaveform.h
 * @brief Independent-source time-domain waveforms (DC, PULSE, SIN, PWL).
 *
 * Faithful to the SPICE/Xyce definitions used by V and I sources. The
 * waveform is evaluated ONCE per time point (not per Newton iteration):
 * the analysis calls loadB(time) each step and the solver keeps the value
 * frozen during the nonlinear iterations of that step — exactly the
 * Xyce architecture (sources load the "B vector", never F).
 */

#ifndef SOURCE_WAVEFORM_H
#define SOURCE_WAVEFORM_H

#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <stdexcept>

namespace e_XSpice {

// C++17 (strict ISO) does not expose M_PI; define the constants we need.
static constexpr double kPi = 3.14159265358979323846;

enum class WaveformType {
    DC,
    PULSE,
    SIN,
    PWL
};

class SourceWaveform {
public:
    SourceWaveform() : type_(WaveformType::DC), dc_value_(0.0) {}

    // ---------------- factories ----------------
    static SourceWaveform dc(double value) {
        SourceWaveform w;
        w.type_ = WaveformType::DC;
        w.dc_value_ = value;
        return w;
    }

    static SourceWaveform pulse(double v1, double v2, double td,
                                double tr, double tf, double pw, double per) {
        SourceWaveform w;
        w.type_ = WaveformType::PULSE;
        w.p_.v1 = v1; w.p_.v2 = v2; w.p_.td = td;
        w.p_.tr = tr; w.p_.tf = tf; w.p_.pw = pw; w.p_.per = per;
        return w;
    }

    static SourceWaveform sin(double v0, double va, double freq,
                              double td, double theta, double phase_deg) {
        SourceWaveform w;
        w.type_ = WaveformType::SIN;
        w.s_.v0 = v0; w.s_.va = va; w.s_.freq = freq;
        w.s_.td = td;  w.s_.theta = theta;
        w.s_.phase = phase_deg * kPi / 180.0;
        return w;
    }

    static SourceWaveform pwl(std::vector<std::pair<double,double>> points) {
        SourceWaveform w;
        w.type_ = WaveformType::PWL;
        w.pwl_ = std::move(points);
        return w;
    }

    // ---------------- evaluation ----------------

    /**
     * @brief Source value at time t.
     */
    double value(double t) const {
        switch (type_) {
            case WaveformType::DC:   return dc_value_;
            case WaveformType::PULSE: return pulseValue(t);
            case WaveformType::SIN:   return sinValue(t);
            case WaveformType::PWL:   return pwlValue(t);
        }
        return dc_value_;
    }

    /**
     * @brief Time points where the waveform has slope discontinuities.
     *        The transient engine must not step over these (Xyce calls
     *        them breakpoints).
     */
    void breakpoints(std::vector<double>& out, double tstop) const {
        switch (type_) {
            case WaveformType::PULSE:
                for (double b : pulseBreakpoints(tstop)) out.push_back(b);
                break;
            case WaveformType::SIN:
                if (s_.va != 0.0 && s_.freq > 0.0 && s_.theta == 0.0) {
                    // one breakpoint per quarter period is enough in practice
                    const double T = 1.0 / s_.freq;
                    for (double t = s_.td; t < tstop; t += T / 4.0) out.push_back(t);
                }
                break;
            case WaveformType::PWL:
                for (const auto& p : pwl_) out.push_back(p.first);
                break;
            case WaveformType::DC:
                break;
        }
    }

    /**
     * @brief Upper bound on the local truncation-safe step size (SPICE
     *        convention: 0.1·PER for PULSE, 0.1/FREQ for SIN).
     */
    double maxTimeStep(double t) const {
        switch (type_) {
            case WaveformType::PULSE:
                return (t < p_.td) ? 0.1 * p_.td
                                   : (p_.per > 0 ? 0.1 * p_.per : 0.0);
            case WaveformType::SIN:
                return (s_.freq > 0.0 && s_.va != 0.0) ? 0.1 / s_.freq : 0.0;
            default:
                return 0.0;  // 0 = "no constraint"
        }
    }

    WaveformType type() const { return type_; }

    // ---------------- AC (small-signal) specification ----------------
    //
    // SPICE gives every source TWO independent specs:
    //     V1 n+ n- DC <value>          -> the DC operating point
    //     V1 n+ n- AC <mag> [phase]    -> the small-signal phasor for .AC
    // The DC waveform above is the time-domain spec; the phasor lives here,
    // next to it, so the device carries both without extra state (exactly
    // what Xyce stores per-source: a DC/B vector plus the AC magnitude and
    // phase read from the netlist).
    //
    // SPICE semantics implemented here:
    //   - "AC mag" without a time spec implies DC 0 (the source contributes
    //     nothing to the operating point),
    //   - phase is in DEGREES and defaults to 0,
    //   - a source with no AC spec is silent in an AC analysis.

    /** @brief Sets the small-signal phasor (mag, phase in degrees). */
    void setAC(double mag, double phase_deg) {
        ac_magnitude_ = mag;
        ac_phase_rad_ = phase_deg * kPi / 180.0;
        ac_given_ = true;
    }

    bool  hasAC() const { return ac_given_; }
    double acMagnitude() const { return ac_magnitude_; }
    double acPhaseRad() const { return ac_phase_rad_; }

    /** @brief Phasor as (cos, sin) — the real/imag parts the AC RHS needs. */
    void acPhasor(double& re, double& im) const {
        re = ac_magnitude_ * std::cos(ac_phase_rad_);
        im = ac_magnitude_ * std::sin(ac_phase_rad_);
    }

private:
    WaveformType type_;
    double dc_value_ = 0.0;
    double ac_magnitude_ = 1.0;
    double ac_phase_rad_ = 0.0;
    bool   ac_given_ = false;

    struct Pulse { double v1, v2, td, tr, tf, pw, per; } p_{};
    struct Sin   { double v0, va, freq, td, theta, phase; } s_{};
    std::vector<std::pair<double,double>> pwl_;

    static constexpr double kEps = 0.0;  // breakpoint tolerance handled by engine

    double pulseValue(double t) const {
        const double td = p_.td, tr = p_.tr, tf = p_.tf,
                     pw = p_.pw, per = p_.per;
        double time = t - td;
        if (time < 0.0) return p_.v1;
        if (per > 0.0 && time > per) {
            const double base = per * std::floor(time / per);
            time -= base;
        }
        if (tr > 0.0 && time < tr)  return p_.v1 + (p_.v2 - p_.v1) * (time / tr);
        if (time < tr + pw)         return p_.v2;
        if (tf > 0.0 && time < tr + pw + tf)
            return p_.v2 + (p_.v1 - p_.v2) * ((time - tr - pw) / tf);
        return p_.v1;
    }

    std::vector<double> pulseBreakpoints(double tstop) const {
        std::vector<double> bp;
        const double td = p_.td, tr = p_.tr, tf = p_.tf, pw = p_.pw;
        bp.push_back(td);
        bp.push_back(td + tr);
        bp.push_back(td + tr + pw);
        bp.push_back(td + tr + pw + tf);
        if (p_.per > 0.0) {
            const double base = p_.per * std::floor((tstop + p_.per) / p_.per);
            bp.push_back(base + td);
            bp.push_back(base + td + p_.per);
            bp.push_back(base + td + p_.per + tr);
            bp.push_back(base + td + p_.per + tr + pw);
            bp.push_back(base + td + p_.per + tr + pw + tf);
        }
        return bp;
    }

    double sinValue(double t) const {
        double time = t - s_.td;
        if (time <= 0.0) {
            return s_.v0 + s_.va * std::sin(s_.phase);
        }
        return s_.v0 + s_.va * std::sin(2.0 * kPi * (s_.freq * time) + s_.phase)
                        * std::exp(-time * s_.theta);
    }

    double pwlValue(double t) const {
        if (pwl_.empty()) return 0.0;
        if (t <= pwl_.front().first) return pwl_.front().second;
        if (t >= pwl_.back().first)  return pwl_.back().second;
        for (size_t i = 1; i < pwl_.size(); ++i) {
            if (t <= pwl_[i].first) {
                const double t1 = pwl_[i-1].first, v1 = pwl_[i-1].second;
                const double t2 = pwl_[i].first,   v2 = pwl_[i].second;
                if (t2 <= t1) return v2;
                return v1 + (v2 - v1) * (t - t1) / (t2 - t1);
            }
        }
        return pwl_.back().second;
    }
};

} // namespace e_XSpice

#endif // SOURCE_WAVEFORM_H
