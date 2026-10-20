/**
 * @author Gabriel Henique Silva
 * @file DaeSolution.h
 * @brief Solution-vector view passed to devices during DAE evaluation.
 *
 * The solution vector x has the MNA layout used throughout e-XSpice:
 *
 *   x = [ v_1 ... v_{n-1} | i_b1 ... i_bk ]
 *
 * where the first block are node voltages (ground is not stored, node 0 is
 * the reference) and the second block are the "branch currents" of
 * voltage-defined devices (independent voltage sources and inductors) —
 * exactly the extra unknowns declared by Device::numExtraUnknowns().
 *
 * Devices read their voltages/currents through this object instead of
 * touching the raw vector, which keeps the index bookkeeping in one place
 * (this is what makes the analysis engines agnostic to the device mix).
 */

#ifndef DAE_SOLUTION_H
#define DAE_SOLUTION_H

#include <complex>
#include <vector>
#include <cassert>

namespace e_XSpice {

class Circuit;
class Node;

class DaeSolution {
public:
    DaeSolution(const std::vector<double>& x,
                const Circuit& circuit,
                const std::vector<int>& node_map)
        : x_(x), im_(nullptr), circuit_(circuit), node_map_(node_map) {}

    /**
     * @brief PHASOR view of the solution (AC small-signal analysis).
     *
     * The imaginary part is optional: an empty vector means "purely real",
     * and every phasor accessor then returns a zero imaginary part. The
     * real-only accessors (voltage/branchValue/data) are untouched, so the
     * hot paths of the DC/transient evaluation pay nothing for this.
     */
    DaeSolution(const std::vector<double>& x_re,
                const std::vector<double>& x_im,
                const Circuit& circuit,
                const std::vector<int>& node_map)
        : x_(x_re), im_(&x_im), circuit_(circuit), node_map_(node_map) {}

    /**
     * @brief Voltage of a node (0.0 for ground).
     * @param node Node object
     * @param fallback If non-null, used when x is smaller than expected
     */
    double voltage(const Node* node) const {
        int idx = rowOf(node);
        return idx >= 0 ? x_[static_cast<size_t>(idx)] : 0.0;
    }

    /**
     * @brief Value of an extra (branch) unknown at absolute index idx.
     */
    double branchValue(int idx) const {
        return idx >= 0 ? x_[static_cast<size_t>(idx)] : 0.0;
    }

    /**
     * @brief Complex voltage of a node (phasor view; 0 for ground).
     */
    std::complex<double> voltagePhasor(const Node* node) const {
        const int idx = rowOf(node);
        if (idx < 0) return std::complex<double>(0.0, 0.0);
        return valueAt(idx);
    }

    /**
     * @brief Complex value of an extra (branch) unknown at absolute index idx.
     */
    std::complex<double> branchPhasor(int idx) const {
        if (idx < 0 || idx >= static_cast<int>(x_.size())) {
            return std::complex<double>(0.0, 0.0);
        }
        return valueAt(idx);
    }

    /**
     * @brief Complex value of the unknown at absolute index idx.
     */
    std::complex<double> valueAt(int idx) const {
        if (idx < 0 || idx >= static_cast<int>(x_.size())) {
            return std::complex<double>(0.0, 0.0);
        }
        const double re = x_[static_cast<size_t>(idx)];
        const double im = (!im_ || idx >= static_cast<int>(im_->size()))
                              ? 0.0
                              : (*im_)[static_cast<size_t>(idx)];
        return std::complex<double>(re, im);
    }

    /**
     * @brief True when the solution carries a non-zero imaginary part.
     */
    bool isPhasor() const { return im_ != nullptr; }

    /**
     * @brief Row index of a node inside x (-1 when ground / unmapped).
     */
    int rowOf(const Node* node) const;

    /**
     * @brief Raw access to the vector (for devices that own contiguous state).
     */
    const std::vector<double>& data() const { return x_; }

    /**
     * @brief Raw access to the imaginary part (empty when purely real).
     */
    const std::vector<double>& imaginaryData() const {
        static const std::vector<double> kEmpty;
        return im_ ? *im_ : kEmpty;
    }

private:
    const std::vector<double>& x_;
    const std::vector<double>* im_;   // null => purely real
    const Circuit& circuit_;
    const std::vector<int>& node_map_;
};

} // namespace e_XSpice

#endif // DAE_SOLUTION_H
