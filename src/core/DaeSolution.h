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
        : x_(x), circuit_(circuit), node_map_(node_map) {}

    /**
     * @brief Voltage of a node (0.0 for ground).
     * @param node Node object
     * @param fallback If non-null, used when x is smaller than expected
     */
    double voltage(const Node* node) const {
        int idx = rowOf(node);
        return idx >= 0 ? x_[idx] : 0.0;
    }

    /**
     * @brief Value of an extra (branch) unknown at absolute index idx.
     */
    double branchValue(int idx) const {
        return idx >= 0 ? x_[idx] : 0.0;
    }

    /**
     * @brief Row index of a node inside x (-1 when ground / unmapped).
     */
    int rowOf(const Node* node) const;

    /**
     * @brief Raw access to the vector (for devices that own contiguous state).
     */
    const std::vector<double>& data() const { return x_; }

private:
    const std::vector<double>& x_;
    const Circuit& circuit_;
    const std::vector<int>& node_map_;
};

} // namespace e_XSpice

#endif // DAE_SOLUTION_H
