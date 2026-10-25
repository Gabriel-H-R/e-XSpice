/**
 * @author Gabriel Henrique Silva
 * @file MeasureEngine.h
 * @brief Post-processing engine for the `.MEASURE`/`.MEAS` and `.FOUR`
 *        commands.
 *
 * The engine is deliberately decoupled from the analysis that produced the
 * data: it takes the time series (or the sweep) and the list of measure
 * definitions, and returns the measured values. That is what makes it
 * reusable across TRAN, DC, AC and NOISE, and what makes the
 * re-measure functionality (re-running the measures on an existing output
 * file) a one-liner.
 *
 * Supported measure types (the most common and useful subset of the Xyce
 * specification):
 *
 *   MAX, MIN, PP, AVG, RMS, INTEG, FIND (AT/WHE N), WHEN, TRIG-TARG,
 *   DERIV (AT/WHEN), FREQ, DUTY, ON_TIME, OFF_TIME, FOUR.
 *
 * The engine works on a UNIFORM representation of the data: a vector of
 * (independent variable, value) pairs, one per variable referenced by the
 * measures. The transient path builds one series per node/branch current;
 * the DC/AC/NOISE paths build one series per sweep point.
 */

#ifndef MEASURE_ENGINE_H
#define MEASURE_ENGINE_H

#include "../core/NetlistParser.h"

#include <functional>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace e_XSpice {

#include "../core/NetlistParser.h"

struct MeasureResult {
    std::string name;
    double value = 0.0;
    double time = 0.0;          // the time/frequency at which it was measured
    bool failed = false;
    std::string message;
};

/**
 * @brief A (time, value) series of a single measured variable.
 */
struct MeasureSeries {
    std::string variable;       // "V(out)", "I(R1)", ...
    std::vector<double> x;      // independent variable (time, freq, sweep)
    std::vector<double> y;      // the variable's value
};

/**
 * @brief Evaluates a `.MEASURE` definition against a series.
 *
 * @param spec    the parsed measure definition
 * @param series  the series of the measured variable (and any TRIG/TARG
 *                variables, looked up by name)
 * @param lookup  a function that returns the series of a variable name
 * @return the measured value
 */
class MeasureEngine {
public:
    static std::vector<MeasureResult> evaluate(
        const std::vector<AnalysisParameters::MeasureSpec>& measures,
        const std::vector<MeasureSeries>& series,
        const std::function<MeasureSeries(const std::string&)>& lookup);

    /**
     * @brief Fourier analysis of a series: DC, fundamental, harmonics, THD.
     *
     * @param freq    the fundamental frequency
     * @param series  the series of the measured variable
     * @param numfreq number of harmonics (including the DC: NUMFREQ)
     * @return a string with the report (the same shape as Xyce's .FOUR
     *         output)
     */
    static std::string fourier(double freq, const MeasureSeries& series,
                                int numfreq = 10);
};

} // namespace e_XSpice

#endif // MEASURE_ENGINE_H
