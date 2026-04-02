/**
 * @author Gabriel Henrique Silva
 * @file CurrentSource.cpp
 * @brief Independent current source implementation
 */

#include "CurrentSource.h"
#include "../solver/SparseMatrix.h"

namespace e_XSpice {

CurrentSource::CurrentSource(const std::string& name,
                             Node* node_plus,
                             Node* node_minus,
                             double current)
    : Device(name, DeviceType::CURRENT_SOURCE)
    , current_(current)
{
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    
    if (node_plus) node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void CurrentSource::stamp(SparseMatrix& G, std::vector<double>& b,
                          const std::vector<int>& node_map) {
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    
    // Current source contributes only to vector b
    // Convention: current flows from n_minus to n_plus
    // 
    // Kirchhoff's law:
    // - Current leaves n_minus: b[j] -= I
    // - Current enters n_plus: b[i] += I
    
    if (i >= 0) {
        b[i] += current_;  // Current enters the positive node
    }
    
    if (j >= 0) {
        b[j] -= current_;  // Current leaves the negative node
    }
    
    // Does not modify matrix G (current source does not affect conductances)
}

double CurrentSource::getPower() const {
    // Power supplied by the source
    // P = (V_plus - V_minus) × I
    double v_plus = nodes_[0]->getVoltage();
    double v_minus = nodes_[1]->getVoltage();
    return (v_plus - v_minus) * current_;
}

} // namespace e_XSpice
