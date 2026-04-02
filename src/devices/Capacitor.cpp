/**
 * @author Gabriel Henrique Silva
 * @file Capacitor.cpp
 * @brief Capacitor implementation
 */

#include "Capacitor.h"
#include "../solver/SparseMatrix.h"
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

Capacitor::Capacitor(const std::string& name,
                     Node* node_plus,
                     Node* node_minus,
                     double capacitance,
                     double initial_voltage)
    : Device(name, DeviceType::CAPACITOR)
    , capacitance_(capacitance)
    , initial_voltage_(initial_voltage)
    , timestep_(0.0)
    , prev_voltage_(initial_voltage)
    , prev_current_(0.0)
{
    if (capacitance <= 0) {
        throw std::invalid_argument("Capacitance must be positive");
    }
    
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    
    if (node_plus) node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void Capacitor::stamp(SparseMatrix& G, std::vector<double>& b,
                     const std::vector<int>& node_map) {
    // For DC analysis: capacitor is open circuit
    // Does not contribute to G or b (DC current = 0)
    (void)G;
    (void)b;
    (void)node_map;
}

void Capacitor::stampTransient(SparseMatrix& G, std::vector<double>& b,
                               const std::vector<int>& node_map,
                               double timestep) {
    timestep_ = timestep;
    
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    
    // Companion model: C is replaced by Geq in parallel with Ieq
    // Geq = 2C/Δt (trapezoidal method) - Disabled for now due to stability issues
    // Ieq = -Geq × V(t-Δt) - I(t-Δt)
    // Geq = C/Δt (backward Euler method)
    // Ieq = Geq × V(t-Δt) (backward Euler, without the I_prev term, for better stability)


    double geq = getEquivalentConductance();
    double ieq = getEquivalentCurrent();
    
    // Stamp equivalent conductance (like a resistor)
    if (i >= 0) {
        G.add(i, i, geq);
        if (j >= 0) {
            G.add(i, j, -geq);
        }
        b[i] += ieq;  // Equivalent current source
    }
    
    if (j >= 0) {
        G.add(j, j, geq);
        if (i >= 0) {
            G.add(j, i, -geq);
        }
        b[j] -= ieq;
    }
}

double Capacitor::getVoltage() const {
    return nodes_[0]->getVoltage() - nodes_[1]->getVoltage();
}

double Capacitor::getCurrent() const {
    if (timestep_ <= 0) {
        return 0.0; // DC
    }
    // Backward Euler: I(t) = Geq × (V(t) - V_prev)
    return getEquivalentConductance() * (getVoltage() - prev_voltage_);
}

/*
double Capacitor::getCurrent() const {
    if (timestep_ <= 0) {
        return 0.0;  // DC
    }
    
    // I = C × dV/dt ≈ C × (V - V_prev) / Δt
    double v = getVoltage();
    return capacitance_ * (v - prev_voltage_) / timestep_;
}
*/

double Capacitor::getPower() const {
    return getVoltage() * getCurrent();
}

double Capacitor::getEnergy() const {
    double v = getVoltage();
    return 0.5 * capacitance_ * v * v;
}

void Capacitor::updateState() {
    // Saves current state for next time step
    prev_voltage_ = getVoltage();
    prev_current_ = getCurrent();
}

double Capacitor::getEquivalentConductance() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Backward Euler: Geq = C / Δt (without the 2 multiplier of Trapezoidal)
    return capacitance_ / timestep_;
}

/*
double Capacitor::getEquivalentConductance() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Trapezoidal method: Geq = 2C/Δt
    return 2.0 * capacitance_ / timestep_;
}
*/

double Capacitor::getEquivalentCurrent() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Backward Euler: Ieq = Geq × V_prev
    // Note that I_prev disappeared, immunizing the method against sudden jumps!
    return getEquivalentConductance() * prev_voltage_;
}

/*
double Capacitor::getEquivalentCurrent() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Ieq = -Geq × V(t-Δt) - I(t-Δt)
    double geq = getEquivalentConductance();
    return -geq * prev_voltage_ - prev_current_;
}
*/

} // namespace e_XSpice
