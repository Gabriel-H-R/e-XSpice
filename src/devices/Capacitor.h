/**
 * @author Gabriel Henrique Silva
 * @file Capacitor.h
 * @brief Capacitor model
 * 
 * Constitutive relation: I = C × dV/dt
 * 
 * For DC analysis: Capacitor is open circuit (I = 0)
 * For AC analysis: Impedance Z = 1/(jωC)
 * For Transient analysis: Uses companion model (temporal discretization)
 * 
 * Companion Model (Trapezoidal Integration Method):
 * At each time step, the capacitor is replaced by:
 * - Equivalent conductance: Geq = 2C/Δt
 * - Equivalent current source: Ieq = -Geq × V(t-Δt) - I(t-Δt)
 */

#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class Capacitor : public Device {
public:
    /**
     * @brief Constructor
     * @param name Capacitor name
     * @param node_plus Positive node
     * @param node_minus Negative node
     * @param capacitance Capacitance (F)
     * @param initial_voltage Initial voltage (V) for transient analysis
     */
    Capacitor(const std::string& name,
              Node* node_plus,
              Node* node_minus,
              double capacitance,
              double initial_voltage = 0.0);
    
    ~Capacitor() override = default;
    
    /**
     * @brief Stamping for DC analysis (open circuit)
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    /**
     * @brief Stamping for transient analysis (companion model)
     * @param G Matrix
     * @param b Vector
     * @param node_map Node mapping
     * @param timestep Time step (Δt)
     */
    void stampTransient(SparseMatrix& G, std::vector<double>& b,
                       const std::vector<int>& node_map,
                       double timestep);
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Calculates current through the capacitor
     * I = C × dV/dt
     */
    double getCurrent() const override;
    
    /**
     * @brief Instantaneous power P = V × I
     */
    double getPower() const override;
    
    /**
     * @brief Stored energy E = 0.5 × C × V²
     */
    double getEnergy() const;
    
    /**
     * @brief Updates state for next time step
     */
    void updateState() override;
    
    /**
     * @brief Sets time step for companion model
     */
    void setTimeStep(double dt) { timestep_ = dt; }
    
    // Getters
    double getCapacitance() const { return capacitance_; }
    double getVoltage() const;
    double getInitialVoltage() const { return initial_voltage_; }

private:
    double capacitance_;         // Capacitance (F)
    double initial_voltage_;     // Initial voltage (V)
    double timestep_;           // Current time step (s)
    
    // Previous state (for integration)
    double prev_voltage_;       // V(t-Δt)
    double prev_current_;       // I(t-Δt)
    
    /**
     * @brief Calculates equivalent conductance
     */
    double getEquivalentConductance() const;
    
    /**
     * @brief Calculates equivalent current (companion model)
     */
    double getEquivalentCurrent() const;
};

} // namespace e_XSpice

#endif // CAPACITOR_H
