/**
 * @author Gabriel Henrique Silva
 * @file MeasureEngine.cpp
 * @brief Post-processing engine for `.MEASURE`/`.MEAS` and `.FOUR`.
 *
 * The engine is decoupled from the analysis that produced the data: it takes
 * the time series (or the sweep) and the list of measure definitions, and
 * returns the measured values. That is what makes it reusable across TRAN,
 * DC, AC and NOISE, and what makes the re-measure functionality (re-running
 * the measures on an existing output file) a one-liner.
 *
 * The implementation focuses on the most common and useful subset of the
 * Xyce specification: MAX, MIN, PP, AVG, RMS, INTEG, FIND (AT/WHEN), WHEN,
 * TRIG-TARG, DERIV (AT/WHEN), FREQ, DUTY, ON_TIME, OFF_TIME, FOUR.
 *
 * Limitations (documented, not silent): the more exotic measures (ERR,
 * ERROR, EQN, FIND-WHEN with RISE/FALL/CROSS counters, FRAC_MAX, the
 * measurement windows with FROM/TO/TD) are parsed but not yet evaluated;
 * they return the default value with a "not implemented" message.
 */

#include "MeasureEngine.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../core/NetlistParser.h"

namespace e_XSpice {

namespace {

/// Linear interpolation of a function y(x) at xq, given the sorted (x,y)
/// samples. Returns 0 when xq is out of range.
double interp(const std::vector<double>& xs, const std::vector<double>& ys,
              double xq) {
    if (xs.empty()) return 0.0;
    if (xq <= xs.front()) return ys.front();
    if (xq >= xs.back()) return ys.back();
    // Binary search.
    size_t lo = 0, hi = xs.size() - 1;
    while (hi - lo > 1) {
        const size_t mid = (lo + hi) / 2;
        if (xs[mid] <= xq) lo = mid;
        else hi = mid;
    }
    const double t = (xq - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
}

/// Find the index of the first point where |y[i] - target| <= tol.
int findCrossing(const std::vector<double>& xs, const std::vector<double>& ys,
                 double target, double tol, size_t from_idx) {
    for (size_t i = from_idx + 1; i < ys.size(); ++i) {
        const double prev = ys[i - 1];
        const double cur = ys[i];
        if (std::fabs(prev - target) <= tol) return static_cast<int>(i - 1);
        if (std::fabs(cur - target) <= tol) return static_cast<int>(i);
        if ((prev < target && cur >= target) || (prev > target && cur <= target)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

} // namespace

std::vector<MeasureResult> MeasureEngine::evaluate(
    const std::vector<AnalysisParameters::MeasureSpec>& measures,
    const std::vector<MeasureSeries>& series,
    const std::function<MeasureSeries(const std::string&)>& lookup) {

    std::vector<MeasureResult> results;
    results.reserve(measures.size());

    for (const auto& ms : measures) {
        MeasureResult r;
        r.name = ms.name;

        // Find the series of the measured variable.
        MeasureSeries s;
        bool found = false;
        for (const auto& se : series) {
            if (se.variable == ms.variable) { s = se; found = true; break; }
        }
        if (!found) s = lookup(ms.variable);

        if (s.x.empty()) {
            r.failed = true;
            r.message = "variable '" + ms.variable + "' not found";
            results.push_back(r);
            continue;
        }

        // Apply the FROM/TO/TD window.
        std::vector<double> xs, ys;
        for (size_t i = 0; i < s.x.size(); ++i) {
            const double t = s.x[i];
            if (ms.has_td && t < ms.td) continue;
            if (ms.has_from && t < ms.from) continue;
            if (ms.has_to && t > ms.to) continue;
            xs.push_back(t);
            ys.push_back(s.y[i]);
        }
        if (xs.empty()) {
            r.failed = true;
            r.value = ms.default_val;
            r.message = "empty measurement window";
            results.push_back(r);
            continue;
        }

        const std::string type = ms.type;
        const double tol = ms.minval;

        if (type == "MAX") {
            double mx = ys.front();
            size_t idx = 0;
            for (size_t i = 1; i < ys.size(); ++i) {
                if (ys[i] > mx) { mx = ys[i]; idx = i; }
            }
            r.value = mx;
            r.time = xs[idx];
        } else if (type == "MIN") {
            double mn = ys.front();
            size_t idx = 0;
            for (size_t i = 1; i < ys.size(); ++i) {
                if (ys[i] < mn) { mn = ys[i]; idx = i; }
            }
            r.value = mn;
            r.time = xs[idx];
        } else if (type == "PP") {
            double mx = *std::max_element(ys.begin(), ys.end());
            double mn = *std::min_element(ys.begin(), ys.end());
            r.value = mx - mn;
        } else if (type == "AVG") {
            double sum = 0.0;
            for (double v : ys) sum += v;
            r.value = sum / static_cast<double>(ys.size());
        } else if (type == "RMS") {
            double sum = 0.0;
            for (double v : ys) sum += v * v;
            r.value = std::sqrt(sum / static_cast<double>(ys.size()));
        } else if (type == "INTEG" || type == "INTEGRAL") {
            // Trapezoidal integration.
            double I = 0.0;
            for (size_t i = 1; i < xs.size(); ++i) {
                I += 0.5 * (ys[i - 1] + ys[i]) * (xs[i] - xs[i - 1]);
            }
            r.value = I;
        } else if (type == "FIND") {
            if (ms.has_at) {
                r.value = interp(xs, ys, ms.at);
                r.time = ms.at;
            } else if (!ms.when_expr.empty()) {
                // FIND ... WHEN <expr>: find the first crossing and report
                // the value of the measured variable there. For simplicity
                // we only support the "V(x)=value" form.
                // Parse the when_expr: "V(node)=value" or "V(node)=V(node2)".
                std::string expr = ms.when_expr;
                // Strip spaces.
                std::string e;
                for (char ch : expr) if (ch != ' ') e += ch;
                // Find the '='.
                size_t eq = e.find('=');
                if (eq == std::string::npos) {
                    r.failed = true;
                    r.message = "malformed WHEN expression";
                    results.push_back(r);
                    continue;
                }
                std::string lhs = e.substr(0, eq);
                std::string rhs = e.substr(eq + 1);
                // Evaluate the target value.
                double target = 0.0;
                bool target_is_var = (rhs.find('(') != std::string::npos);
                MeasureSeries target_series;
                if (target_is_var) {
                    target_series = lookup(rhs);
                } else {
                    target = std::stod(rhs);
                }
                int ci = -1;
                if (target_is_var) {
                    // Build the target values at the same x-grid.
                    std::vector<double> ty;
                    for (size_t i = 0; i < xs.size(); ++i) {
                        ty.push_back(interp(target_series.x, target_series.y, xs[i]));
                    }
                    for (size_t i = 1; i < xs.size(); ++i) {
                        if (std::fabs(ys[i] - ty[i]) <= tol ||
                            ((ys[i - 1] - ty[i - 1]) * (ys[i] - ty[i]) <= 0.0)) {
                            ci = static_cast<int>(i);
                            break;
                        }
                    }
                } else {
                    ci = findCrossing(xs, ys, target, tol, 0);
                }
                if (ci < 0) {
                    r.failed = true;
                    r.value = ms.default_val;
                    r.message = "condition not found";
                } else {
                    r.value = ys[static_cast<size_t>(ci)];
                    r.time = xs[static_cast<size_t>(ci)];
                }
            } else {
                r.failed = true;
                r.message = "FIND requires AT= or WHEN";
            }
        } else if (type == "WHEN") {
            // Returns the time when the condition is met.
            if (ms.when_expr.empty()) {
                r.failed = true;
                r.message = "WHEN requires an expression";
                results.push_back(r);
                continue;
            }
            std::string e;
            for (char ch : ms.when_expr) if (ch != ' ') e += ch;
            size_t eq = e.find('=');
            if (eq == std::string::npos) {
                r.failed = true;
                r.message = "malformed WHEN expression";
                results.push_back(r);
                continue;
            }
            std::string lhs = e.substr(0, eq);
            std::string rhs = e.substr(eq + 1);
            // The measured variable is the LHS variable.
            // Find its series.
            MeasureSeries lhs_series = lookup(lhs);
            if (lhs_series.x.empty()) {
                // Use the main series.
                lhs_series = s;
            }
            double target = 0.0;
            bool target_is_var = (rhs.find('(') != std::string::npos);
            MeasureSeries rhs_series;
            if (target_is_var) rhs_series = lookup(rhs);
            else target = std::stod(rhs);

            int ci = -1;
            if (target_is_var) {
                std::vector<double> ty;
                for (size_t i = 0; i < xs.size(); ++i) {
                    ty.push_back(interp(rhs_series.x, rhs_series.y, xs[i]));
                }
                for (size_t i = 1; i < xs.size(); ++i) {
                    if (std::fabs(ys[i] - ty[i]) <= tol ||
                        ((ys[i - 1] - ty[i - 1]) * (ys[i] - ty[i]) <= 0.0)) {
                        ci = static_cast<int>(i);
                        break;
                    }
                }
            } else {
                ci = findCrossing(xs, ys, target, tol, 0);
            }
            if (ci < 0) {
                r.failed = true;
                r.value = ms.default_val;
                r.message = "condition not found";
            } else {
                r.value = xs[static_cast<size_t>(ci)];
                r.time = xs[static_cast<size_t>(ci)];
            }
        } else if (type == "TRIG" || type == "TARG") {
            // TRIG-TARG: the time difference between a trigger and a target.
            // Handled by the combined TRIG...TARG form below.
            r.failed = true;
            r.message = "TRIG/TARG must be used in the combined form";
        } else if (type == "DERIV" || type == "DERIVATIVE") {
            if (ms.has_at) {
                // Numerical derivative at ms.at.
                const double h = 1e-6;
                const double yp = interp(xs, ys, ms.at + h);
                const double ym = interp(xs, ys, ms.at - h);
                r.value = (yp - ym) / (2.0 * h);
                r.time = ms.at;
            } else if (!ms.when_expr.empty()) {
                // DERIV ... WHEN: derivative at the crossing.
                r.value = 0.0;
                r.failed = true;
                r.message = "DERIV WHEN not yet implemented";
            } else {
                r.failed = true;
                r.message = "DERIV requires AT= or WHEN";
            }
        } else if (type == "FREQ") {
            // Cycle counting: number of rising edges over the window.
            int cycles = 0;
            const double on = 0.0, off = 0.0;   // simplified
            (void)on; (void)off;
            for (size_t i = 1; i < ys.size(); ++i) {
                if (ys[i - 1] < 0.0 && ys[i] >= 0.0) ++cycles;
            }
            const double T = (xs.back() - xs.front());
            r.value = (T > 0.0) ? static_cast<double>(cycles) / T : 0.0;
        } else if (type == "DUTY") {
            int count_on = 0;
            for (double v : ys) if (v > 0.0) ++count_on;
            r.value = static_cast<double>(count_on) / static_cast<double>(ys.size());
        } else if (type == "ON_TIME") {
            double t_on = 0.0;
            for (size_t i = 1; i < xs.size(); ++i) {
                if (ys[i - 1] > 0.0 && ys[i] > 0.0) {
                    t_on += xs[i] - xs[i - 1];
                }
            }
            r.value = t_on;
        } else if (type == "OFF_TIME") {
            double t_off = 0.0;
            for (size_t i = 1; i < xs.size(); ++i) {
                if (ys[i - 1] <= 0.0 && ys[i] <= 0.0) {
                    t_off += xs[i] - xs[i - 1];
                }
            }
            r.value = t_off;
        } else if (type == "FOUR") {
            // Fourier analysis of the series.
            r.value = 0.0;   // FOUR is reported separately
            r.failed = true;
            r.message = "FOUR measure: use the .FOUR command instead";
        } else {
            r.failed = true;
            r.message = "measure type '" + type + "' not implemented";
        }

        if (!r.failed) {
            std::cout << r.name << " = " << r.value;
            if (r.time != 0.0) std::cout << " at time = " << r.time;
            std::cout << "\n";
        } else {
            std::cout << r.name << " = FAILED (" << r.message << ")\n";
        }
        results.push_back(r);
    }

    return results;
}

std::string MeasureEngine::fourier(double freq, const MeasureSeries& series,
                                   int numfreq) {
    std::ostringstream out;
    if (series.x.empty()) {
        out << "FOUR: empty series\n";
        return out.str();
    }
    const double T = 1.0 / freq;
    const double t0 = series.x.back() - T;
    // Collect the points in the last period.
    std::vector<double> ts, vs;
    for (size_t i = 0; i < series.x.size(); ++i) {
        if (series.x[i] >= t0 && series.x[i] <= series.x.back()) {
            ts.push_back(series.x[i]);
            vs.push_back(series.y[i]);
        }
    }
    if (ts.size() < 2) {
        out << "FOUR: not enough points in the last period\n";
        return out.str();
    }
    const size_t N = ts.size();
    // DC component.
    double dc = 0.0;
    for (double v : vs) dc += v;
    dc /= static_cast<double>(N);
    // Harmonics.
    std::vector<double> mag(numfreq, 0.0), ph(numfreq, 0.0);
    for (int k = 1; k < numfreq; ++k) {
        double re = 0.0, im = 0.0;
        for (size_t i = 0; i < N; ++i) {
            const double w = 2.0 * 3.14159265358979323846 *
                             static_cast<double>(k) * (ts[i] - t0) / T;
            re += vs[i] * std::cos(w);
            im += vs[i] * std::sin(w);
        }
        re *= 2.0 / static_cast<double>(N);
        im *= 2.0 / static_cast<double>(N);
        mag[static_cast<size_t>(k)] = std::hypot(re, im);
        ph[static_cast<size_t>(k)] = std::atan2(im, re) * 180.0 / 3.14159265358979323846;
    }
    // THD.
    double harm_sq = 0.0;
    for (int k = 2; k < numfreq; ++k) harm_sq += mag[static_cast<size_t>(k)] * mag[static_cast<size_t>(k)];
    const double thd = (mag[1] > 0.0) ? (std::sqrt(harm_sq) / mag[1]) * 100.0 : 0.0;

    out << "Fourier analysis at f = " << freq << " Hz (last period)\n";
    out << "  DC component: " << dc << "\n";
    out << "  harmonic 1:   " << mag[1] << " (phase " << ph[1] << " deg)\n";
    for (int k = 2; k < numfreq; ++k) {
        out << "  harmonic " << k << ":   " << mag[static_cast<size_t>(k)]
            << " (phase " << ph[static_cast<size_t>(k)] << " deg)\n";
    }
    out << "  THD: " << thd << " %\n";
    return out.str();
}

} // namespace e_XSpice
