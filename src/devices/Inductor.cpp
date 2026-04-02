/**
 * @author Gabriel Henrique Silva
 * @file Inductor.cpp
 * @brief Inductor implementation
 */

#include "Inductor.h"
#include "../solver/SparseMatrix.h"
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

Inductor::Inductor(const std::string& name,
                   Node* node_plus,
                   Node* node_minus,
                   double inductance,
                   double initial_current)
    : Device(name, DeviceType::INDUCTOR)
    , inductance_(inductance)
    , initial_current_(initial_current)
    , timestep_(0.0)
    , current_(initial_current)
    , current_index_(-1)
    , prev_voltage_(0.0)
    , prev_current_(initial_current)
{
    if (inductance <= 0) {
        throw std::invalid_argument("Inductance must be positive");
    }
    
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    
    if (node_plus) node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void Inductor::stamp(SparseMatrix& G, std::vector<double>& b,
                    const std::vector<int>& node_map) {
    // For DC analysis: inductor is short circuit (V = 0)
    // Similar to voltage source with V = 0
    
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    int k = current_index_;
    
    if (k < 0) {
        throw std::runtime_error("Current index not set for inductor " + name_);
    }
    
    // Stamp as voltage source with V = 0
    if (i >= 0) {
        G.add(i, k, 1.0);
        G.add(k, i, 1.0);
    }
    
    if (j >= 0) {
        G.add(j, k, -1.0);
        G.add(k, j, -1.0);
    }
    
    b[k] = 0.0;  // V = 0 (DC short circuit)
}

void Inductor::stampTransient(SparseMatrix& G, std::vector<double>& b,
                              const std::vector<int>& node_map,
                              double timestep) {
    timestep_ = timestep;
    
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    
    // Norton equivalent companion model:
    // Geq = Δt / (2L)
    // Ieq = I(t-Δt) + Geq × V(t-Δt)
    
    double geq = getEquivalentConductance();
    double ieq = getEquivalentCurrent();
    
    // Stamp conductance
    if (i >= 0) {
        G.add(i, i, geq);
        if (j >= 0) {
            G.add(i, j, -geq);
        }
        b[i] += ieq;
    }
    
    if (j >= 0) {
        G.add(j, j, geq);
        if (i >= 0) {
            G.add(j, i, -geq);
        }
        b[j] -= ieq;
    }
}

double Inductor::getVoltage() const {
    return nodes_[0]->getVoltage() - nodes_[1]->getVoltage();
}

double Inductor::getPower() const {
    return getVoltage() * current_;
}

double Inductor::getEnergy() const {
    return 0.5 * inductance_ * current_ * current_;
}

void Inductor::updateState() {
    // For transient: calculates new current
    if (timestep_ > 0) {
        // I = I(t-Δt) + (1/L) × integral[V(τ)dτ]
        // Using trapezoidal: I = I_prev + (Δt/2L) × (V + V_prev)
        double v = getVoltage();
        current_ = prev_current_ + (timestep_ / (2.0 * inductance_)) * (v + prev_voltage_);
    }
    
    // Save state
    prev_voltage_ = getVoltage();
    prev_current_ = current_;
}

double Inductor::getEquivalentConductance() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Norton: Geq = Δt / (2L)
    return timestep_ / (2.0 * inductance_);
}

double Inductor::getEquivalentCurrent() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Norton: Ieq = I(t-Δt) + Geq × V(t-Δt)
    double geq = getEquivalentConductance();
    return prev_current_ + geq * prev_voltage_;
}

} // namespace e_XSpice
