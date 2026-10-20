/**
 * @author Gabriel Henrique Silva
 * @file PZAnalysis.h
 * @brief Pole-zero analysis (`.PZ`): poles and zeros of the small-signal
 *        transfer function, from a single linearization at the operating point.
 *
 * Card (ngspice / SPICE3)
 * -----------------------
 *     .PZ <node1> <node2> <node3> <node4> (vol|cur) (pol|zer|pz)
 *
 *   - node1, node2: output nodes (V(node1) − V(node2))
 *   - node3, node4: input nodes  (V(node3) − V(node4) for `vol`,
 *     or the current through the source between them for `cur`)
 *   - `vol`: transfer function = V(out)/V(in)
 *   - `cur`: transfer function = V(out)/I(in)
 *   - `pol`: poles only, `zer`: zeros only, `pz`: both
 *
 * Method
 * ------
 * At the DC operating point the small-signal system is `A(s) = G + s·C`.
 * The transfer function is `H(s) = cᵀ · A(s)⁻¹ · b`.
 *
 * - Poles are the finite eigenvalues of the generalized pencil `(G, C)`:
 *   `G·v = λ·C·v`. They are found in ONE shot by forming `M = G⁻¹·C`
 *   (N solves of G on the existing sparse LU) and running a QR algorithm
 *   with Wilkinson shifts on M. Cost: O(N³) for all poles, versus the
 *   O(N⁴) iterative search that SPICE3/ngspice use and that is known to
 *   "get lost" on large circuits.
 * - Zeros are found by Newton's method on `H(s)` and `H'(s)`, where both
 *   are computed from two solves per iteration (forward + adjoint) — the
 *   same trick as the noise analysis. Multiple starting points cover the
 *   complex plane; deflation removes found zeros.
 * - The DC gain `H(0) = cᵀ·G⁻¹·b` is reported too.
 *
 * This is the standard approach in modern (MATLAB/Octave) control toolboxes,
 * adapted to the sparse DAE infrastructure: the factorization is reused for
 * every solve, and the matrix is the same one the AC and noise analyses use.
 *
 * Note: Xyce does NOT implement `.PZ` at all (its reference guide lists it
 * as unsupported). This implementation is based on ngspice/SPICE3 semantics.
 */

#ifndef PZ_ANALYSIS_H
#define PZ_ANALYSIS_H

#include "../core/Circuit.h"

#include <complex>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief One pole or zero of the transfer function.
 */
struct PZPoint {
    double real = 0.0;
    double imag = 0.0;
    bool is_pole = true;
};

struct PZResult {
    bool success = false;
    std::string message;

    std::string output_nodes;   // "node1,node2"
    std::string input_nodes;    // "node3,node4"
    bool current_input = false; // true = V(out)/I(in), false = V(out)/V(in)

    double dc_gain = 0.0;       // H(0) = cᵀ·G⁻¹·b

    std::vector<PZPoint> poles;
    std::vector<PZPoint> zeros;
};

class PZAnalysis {
public:
    explicit PZAnalysis(Circuit& circuit);
    ~PZAnalysis();

    void setOutputNodes(const std::string& n1, const std::string& n2) {
        out1_ = n1; out2_ = n2;
    }
    void setInputNodes(const std::string& n3, const std::string& n4) {
        in1_ = n3; in2_ = n4;
    }
    void setCurrentInput(bool on) { current_input_ = on; }
    void setPolesOnly(bool on) { poles_only_ = on; }
    void setZerosOnly(bool on) { zeros_only_ = on; }

    void setMaxIterations(int n) { max_iterations_ = n; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setVerbose(bool v) { verbose_ = v; }

    PZResult run();

private:
    Circuit& circuit_;
    std::string out1_, out2_, in1_, in2_;
    bool current_input_ = false;
    bool poles_only_ = false;
    bool zeros_only_ = false;

    int    max_iterations_ = 200;
    double tolerance_ = 1e-9;
    bool   verbose_ = false;
};

} // namespace e_XSpice

#endif // PZ_ANALYSIS_H