/**
 * @file Resistor.cpp
 * @brief Implementação do resistor linear
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
    
    // Conecta nós
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    
    // Registra este dispositivo nos nós
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
    // Obtém nós
    Node* n_plus = nodes_[0];
    Node* n_minus = nodes_[1];
    
    // Obtém índices na matriz (ground é -1)
    int i = n_plus->isGround() ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    
    double g = conductance_;
    
    // Stamp da condutância na matriz G
    // Resistor entre nós i e j contribui:
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
    
    // Resistor não contribui para o vetor b (sem fontes independentes)
}

void Resistor::stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                       std::vector<std::complex<double>>& b_ac,
                       const std::vector<int>& node_map,
                       double omega) {
    (void)omega; // Resistor não depende de frequência
    (void)b_ac;
    Node* n_plus  = nodes_[0];
    Node* n_minus = nodes_[1];
    int i = n_plus->isGround()  ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    std::complex<double> g(conductance_, 0.0); // Y = G (real)
    if (i >= 0) { Z[i][i] += g; if (j >= 0) Z[i][j] -= g; }
    if (j >= 0) { Z[j][j] += g; if (i >= 0) Z[j][i] -= g; }
}

double Resistor::getCurrent() const {
    // Corrente de n_plus para n_minus
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