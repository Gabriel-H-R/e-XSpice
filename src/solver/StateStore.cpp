/**
 * @author Gabriel Henrique Silva
 * @file StateStore.cpp
 * @brief Owned-by-device solution/state history store implementation.
 */

#include "StateStore.h"

namespace e_XSpice {

StateIndex StateStore::registerStates(int count) {
    StateIndex idx;
    idx.base  = total_size_;
    idx.count = count;
    total_size_ += count;

    // Grow all history vectors so they stay the same length.
    q_curr_.resize(total_size_, 0.0);
    q_next_.resize(total_size_, 0.0);
    qdot_curr_.resize(total_size_, 0.0);

    return idx;
}

void StateStore::setDeviceQ(const StateIndex& idx, const double* values, int count) {
    if (count != idx.count) {
        throw std::invalid_argument("StateStore::setDeviceQ: count mismatch");
    }
    std::copy(values, values + count, q_next_.begin() + idx.base);
}

void StateStore::acceptStep(double s, bool trapezoidal) {
    for (int i = 0; i < total_size_; ++i) {
        // dQ/dt at the new point from the integrator formula
        double qdot_new = s * (q_next_[i] - q_curr_[i]);
        if (trapezoidal) {
            qdot_new -= qdot_curr_[i];   // trapezoidal history
        }
        qdot_curr_[i] = qdot_new;        // accepted secondary state
        q_curr_[i]    = q_next_[i];      // rotate primary history
    }
}

} // namespace e_XSpice
