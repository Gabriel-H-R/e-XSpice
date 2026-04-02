/**
 * @author Gabriel Henrique Silva
 * @file Inductor.h
 * @brief Inductor model
 * 
 * Constitutive relation: V = L × dI/dt
 * 
 * For DC analysis: Inductor is short circuit (V = 0)
 * For AC analysis: Impedance Z = jωL
 * For Transient analysis: Uses companion model
 * 
 * Companion Model (Trapezoidal Integration Method):
 * The inductor is replaced by:
 * - Equivalent resistance: Req = 2L/Δt
 * - Equivalent voltage source in series
 * 
 * Since voltage sources add variables, we use Norton equivalent:
 * - Conductance: Geq = Δt/(2L)
 * - Current source: Ieq = I(t-Δt) + Geq × V(t-Δt)
 */

#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class Inductor : public Device {
public:
    /**
     * @brief Constructor
     * @param name Inductor name
     * @param node_plus Positive node
     * @param node_minus Negative node
     * @param inductance Inductance (H)
     * @param initial_current Initial current (A) for transient analysis
     */
    Inductor(const std::string& name,
             Node* node_plus,
             Node* node_minus,
             double inductance,
             double initial_current = 0.0);
    
    ~Inductor() override = default;
    
    /**
     * @brief Stamping for DC analysis
     * DC: inductor is short circuit, but needs current variable
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    /**
     * @brief Stamping for transient analysis (companion model)
     */
    void stampTransient(SparseMatrix& G, std::vector<double>& b,
                       const std::vector<int>& node_map,
                       double timestep);
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Returns current through the inductor
     */
    double getCurrent() const override { return current_; }
    
    /**
     * @brief Instantaneous power P = V × I
     */
    double getPower() const override;
    
    /**
     * @brief Stored energy E = 0.5 × L × I²
     */
    double getEnergy() const;
    
    /**
     * @brief Updates state for next time step
     */
    void updateState() override;
    
    /**
     * @brief Sets time step
     */
    void setTimeStep(double dt) { timestep_ = dt; }
    
    /**
     * @brief Sets current index in the matrix (like voltage source)
     */
    void setCurrentIndex(int idx) { current_index_ = idx; }
    int getCurrentIndex() const { return current_index_; }
    
    // Getters
    double getInductance() const { return inductance_; }
    double getVoltage() const;
    double getInitialCurrent() const { return initial_current_; }

private:
    double inductance_;          // Inductance (H)
    double initial_current_;     // Initial current (A)
    double timestep_;           // Current time step (s)
    double current_;            // Current current (A)
    int current_index_;         // Index in expanded matrix (for DC)
    
    // Previous state (for integration)
    double prev_voltage_;       // V(t-Δt)
    double prev_current_;       // I(t-Δt)
    
    /**
     * @brief Calculates equivalent conductance (Norton)
     */
    double getEquivalentConductance() const;
    
    /**
     * @brief Calculates equivalent current (Norton)
     */
    double getEquivalentCurrent() const;
};

} // namespace e_XSpice

#endif // INDUCTOR_H
