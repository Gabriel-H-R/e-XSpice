/**
 * @author Gabriel Henrique Silva
 * @file StateStore.h
 * @brief Owned-by-device solution/state history store for the DAE formulation.
 *
 * Xyce analogy:
 * - The state of a charge-storing device is split into PRIMARY state (Q:
 *   charges/fluxes, stored per accepted time point) and SECONDARY state
 *   (dQ/dt, derived by the integrator from the primary history).
 * - Each device registers how many states it owns and receives a
 *   contiguous range; indices never change for the simulation lifetime.
 *
 * Improvements over a naive implementation:
 * - Devices OWN their state ranges, so history advance/restore needs no
 *   per-device maps and no per-step allocation.
 * - One flat vector per history level (better cache locality than one
 *   std::vector per device); rotating history is a linear pass, not a copy
 *   of every double per device.
 * - Save/restore of the accepted point (for LTE step rejection) is
 *   O(#states) and never reallocates.
 */

#ifndef STATE_STORE_H
#define STATE_STORE_H

#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include "../core/StateIndex.h"

namespace e_XSpice {

class StateStore {
public:
    StateStore() = default;

    /**
     * @brief Registers `count` states for a device.
     * @return handle the device keeps for the simulation lifetime
     */
    StateIndex registerStates(int count);

    /**
     * @brief Total registered states.
     */
    int size() const { return total_size_; }

    /**
     * @brief Writes the candidate Q values for one device.
     */
    void setDeviceQ(const StateIndex& idx, const double* values, int count);

    /**
     * @brief Accepted Q of a device's own state.
     */
    double deviceQ(const StateIndex& idx, int local = 0) const {
        checkSlot(idx.base + local);
        return q_curr_[idx.base + local];
    }

    /**
     * @brief Accepted dQ/dt of a device's own state (capacitor current).
     */
    double deviceQDot(const StateIndex& idx, int local = 0) const {
        checkSlot(idx.base + local);
        return qdot_curr_[idx.base + local];
    }

    /**
     * @brief Initializes the history: the accepted charge becomes the
     *        current candidate and dQ/dt is reset to zero (start of a
     *        transient, after an operating-point solve).
     */
    void initializeHistory() {
        for (int i = 0; i < total_size_; ++i) {
            q_curr_[i]    = q_next_[i];
            qdot_curr_[i] = 0.0;
        }
    }

    /**
     * @brief Advances every state by one accepted time step.
     *
     * Computes the secondary state (dQ/dt) with the integrator formula and
     * rotates the primary history (Xyce: updatePrimaryState + the
     * integrator's updateStateDeriv):
     *
     *   BE :     dQ/dt = s·(Q_next − Q_curr)                    (s = 1/dt)
     *   TRAP:    dQ/dt = s·(Q_next − Q_curr) − dQ/dt_prev       (s = 2/dt)
     *
     * @param s             integration coefficient
     * @param trapezoidal   trapezoidal reconstruction when true
     */
    void acceptStep(double s, bool trapezoidal);

    /**
     * @brief Snapshot of the accepted point (O(#states) copy).
     */
    void saveAccepted() {
        q_curr_saved_    = q_curr_;
        qdot_curr_saved_ = qdot_curr_;
    }

    /**
     * @brief Restores the last snapshot (after a rejected time step).
     */
    void restoreAccepted() {
        q_curr_    = q_curr_saved_;
        qdot_curr_ = qdot_curr_saved_;
    }

    /**
     * @brief Resets every state (used when a new simulation starts).
     */
    void reset() {
        std::fill(q_curr_.begin(), q_curr_.end(), 0.0);
        std::fill(q_next_.begin(), q_next_.end(), 0.0);
        std::fill(qdot_curr_.begin(), qdot_curr_.end(), 0.0);
    }

private:
    int total_size_ = 0;

    std::vector<double> q_curr_;     // Q(t_n)   (accepted)
    std::vector<double> q_next_;     // Q(t_{n+1}) (candidate)
    std::vector<double> qdot_curr_;  // dQ/dt(t_n) (accepted)

    std::vector<double> q_curr_saved_;   // snapshot (empty until used)
    std::vector<double> qdot_curr_saved_;

    void checkSlot(int slot) const {
        if (slot < 0 || slot >= total_size_) {
            throw std::out_of_range("StateStore: invalid state slot");
        }
    }
};

} // namespace e_XSpice

#endif // STATE_STORE_H
