/**
 * @author Gabriel Henrique Silva
 * @file Resistor.cpp
 * @brief Linear resistor implementation
 */

#include "Resistor.h"
#include "../solver/SparseMatrix.h"
#include <stdexcept>
#include <cmath>

namespace e_XSpice {

Resistor::Resistor(const std::string& name, 
                   Node* node_plus, 
                   Node* node_minus,
                   double resistance)
    : Device(name, DeviceType::RESISTOR)
    , resistance_(resistance)
{
    validateResistance(resistance);
    
    conductance_ = 1.0 / resistance_;
    
    // Connects nodes
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    
    // Registers this device in the nodes
    if (node_plus) node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void Resistor::validateResistance(double r) {
    if (r <= 0.0) {
        throw std::invalid_argument("Resistance must be positive");
    }
    if (std::isinf(r) || std::isnan(r)) {
        throw std::invalid_argument("Resistance must be finite");
    }
}

void Resistor::stamp(SparseMatrix& G, std::vector<double>& b,
                     const std::vector<int>& node_map) {
    // Gets nodes
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    // Gets indices in matrix (ground is -1)
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    
    double g = conductance_;
    
    // Stamping conductance in matrix G
    // Resistor between nodes i and j contributes:
    // G[i][i] += g,  G[i][j] -= g
    // G[j][i] -= g,  G[j][j] += g
    
    if (i >= 0) {
        G.add(i, i, g);  // Diagonal
        if (j >= 0) {
            G.add(i, j, -g);  // Off-diagonal
        }
    }
    
    if (j >= 0) {
        G.add(j, j, g);  // Diagonal
        if (i >= 0) {
            G.add(j, i, -g);  // Off-diagonal
        }
    }
    
    // Resistor does not contribute to vector b (no independent sources)
}

double Resistor::getCurrent() const {
    // Current from n_plus to n_minus
    // I = (V_plus - V_minus) / R = (V_plus - V_minus) * g
    double v_plus = nodes_[0]->getVoltage();
    double v_minus = nodes_[1]->getVoltage();
    return (v_plus - v_minus) * conductance_;
}

double Resistor::getPower() const {
    // P = I² × R
    double current = getCurrent();
    return current * current * resistance_;
}

} // namespace e_XSpice
