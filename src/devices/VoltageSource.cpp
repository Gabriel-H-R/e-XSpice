/**
 * @author Gabriel Henrique Silva
 * @file VoltageSource.cpp
 * @brief Independent voltage source implementation
 */

#include "VoltageSource.h"
#include "../solver/SparseMatrix.h"
#include <stdexcept>

namespace e_XSpice {

VoltageSource::VoltageSource(const std::string& name,
                             Node* node_plus,
                             Node* node_minus,
                             double voltage)
    : Device(name, DeviceType::VOLTAGE_SOURCE)
    , voltage_(voltage)
    , current_(0.0)
    , current_index_(-1)
{
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    
    if (node_plus) node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void VoltageSource::stamp(SparseMatrix& G, std::vector<double>& b,
                          const std::vector<int>& node_map) {
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    int k = current_index_;  // Current index in the expanded matrix
    
    if (k < 0) {
        throw std::runtime_error("Current index not set for voltage source " + name_);
    }
    
    // Contribution to matrix G (expanded system)
    // 
    // Node equations (KCL):
    // If current flows from + to -, then:
    // - Current leaves node i: G[i][k] += 1
    // - Current enters node j: G[j][k] -= 1
    //
    // Source voltage equation:
    // V_i - V_j = V_source
    // G[k][i] = 1, G[k][j] = -1, b[k] = V_source
    
    if (i >= 0) {
        G.add(i, k, 1.0);   // Current leaves i
        G.add(k, i, 1.0);   // Voltage of i
    }
    
    if (j >= 0) {
        G.add(j, k, -1.0);  // Current enters j
        G.add(k, j, -1.0);  // Voltage of j
    }

    G.add(k, k, -1e-12); // Regularization term to avoid singularity (Remove later?)
    
    // Right-hand side: source voltage
    b[k] = voltage_;
}

double VoltageSource::getPower() const {
    // Power supplied = V × I
    // If I > 0, the source is supplying energy
    // If I < 0, the source is absorbing energy
    return voltage_ * current_;
}

} // namespace e_XSpice
