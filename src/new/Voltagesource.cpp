/**
 * @file VoltageSource.cpp
 * @brief Implementação da fonte de tensão independente
 */

#include "VoltageSource.h"
#include "../solver/SparseMatrix.h"
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace e_XSpice {

VoltageSource::VoltageSource(const std::string& name,
                             Node* node_plus,
                             Node* node_minus,
                             double voltage)
    : Device(name, DeviceType::VOLTAGE_SOURCE)
    , voltage_(voltage)
    , current_(0.0)
    , current_index_(-1)
    , ac_magnitude_(1.0)   // Padrão: 1V AC (amplitude de análise small-signal)
    , ac_phase_rad_(0.0)
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
    int k = current_index_;  // Índice da corrente na matriz expandida
    
    if (k < 0) {
        throw std::runtime_error("Current index not set for voltage source " + name_);
    }
    
    // Contribuição para a matriz G (sistema expandido)
    // 
    // Equações de nó (KCL):
    // Se a corrente flui de + para -, então:
    // - Corrente sai do nó i: G[i][k] += 1
    // - Corrente entra no nó j: G[j][k] -= 1
    //
    // Equação de tensão da fonte:
    // V_i - V_j = V_source
    // G[k][i] = 1, G[k][j] = -1, b[k] = V_source
    
    if (i >= 0) {
        G.add(i, k, 1.0);   // Corrente sai de i
        G.add(k, i, 1.0);   // Tensão de i
    }
    
    if (j >= 0) {
        G.add(j, k, -1.0);  // Corrente entra em j
        G.add(k, j, -1.0);  // Tensão de j
    }
    
    // Lado direito: tensão da fonte
    b[k] = voltage_;
}

void VoltageSource::stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                            std::vector<std::complex<double>>& b_ac,
                            const std::vector<int>& node_map,
                            double omega) {
    (void)omega;
    Node* n_plus  = nodes_[0];
    Node* n_minus = nodes_[1];
    int i = n_plus->isGround()  ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    int k = current_index_;
    if (k < 0 || k >= static_cast<int>(b_ac.size())) return;

    // Mesma estrutura KVL do DC, mas RHS é amplitude AC complexa
    std::complex<double> V_ac = std::polar(ac_magnitude_, ac_phase_rad_);

    if (i >= 0) { Z[i][k] += std::complex<double>(1.0,0); Z[k][i] += std::complex<double>(1.0,0); }
    if (j >= 0) { Z[j][k] -= std::complex<double>(1.0,0); Z[k][j] -= std::complex<double>(1.0,0); }
    b_ac[k] = V_ac;
}

double VoltageSource::getPower() const {
    // Potência fornecida = V × I
    // Se I > 0, a fonte está fornecendo energia
    // Se I < 0, a fonte está absorvendo energia
    return voltage_ * current_;
}

} // namespace e_XSpice