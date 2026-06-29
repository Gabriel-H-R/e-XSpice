/**
 * @file Capacitor.cpp
 * @brief Implementação do capacitor
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
    // Para análise DC: capacitor é circuito aberto
    // Não contribui para G ou b (corrente DC = 0)
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
    
    // Companion model: C é substituído por Geq em paralelo com Ieq
    // Geq = 2C/Δt (método trapezoidal)
    // Ieq = -Geq × V(t-Δt) - I(t-Δt)
    
    double geq = getEquivalentConductance();
    double ieq = getEquivalentCurrent();
    
    // Stamp condutância equivalente (como resistor)
    if (i >= 0) {
        G.add(i, i, geq);
        if (j >= 0) {
            G.add(i, j, -geq);
        }
        b[i] += ieq;  // Fonte de corrente equivalente
    }
    
    if (j >= 0) {
        G.add(j, j, geq);
        if (i >= 0) {
            G.add(j, i, -geq);
        }
        b[j] -= ieq;
    }
}

void Capacitor::stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                        std::vector<std::complex<double>>& b_ac,
                        const std::vector<int>& node_map,
                        double omega) {
    (void)b_ac;
    Node* n_plus  = nodes_[0];
    Node* n_minus = nodes_[1];
    int i = n_plus->isGround()  ? -1 : node_map[n_plus->getId()];
    int j = n_minus->isGround() ? -1 : node_map[n_minus->getId()];
    // Y_C = jωC — imaginário puro. j em C++ = complex(0,1)
    std::complex<double> y(0.0, omega * capacitance_);
    if (i >= 0) { Z[i][i] += y; if (j >= 0) Z[i][j] -= y; }
    if (j >= 0) { Z[j][j] += y; if (i >= 0) Z[j][i] -= y; }
}

double Capacitor::getVoltage() const {
    return nodes_[0]->getVoltage() - nodes_[1]->getVoltage();
}

double Capacitor::getCurrent() const {
    if (timestep_ <= 0) {
        return 0.0;  // DC
    }
    
    // I = C × dV/dt ≈ C × (V - V_prev) / Δt
    double v = getVoltage();
    return capacitance_ * (v - prev_voltage_) / timestep_;
}

double Capacitor::getPower() const {
    return getVoltage() * getCurrent();
}

double Capacitor::getEnergy() const {
    double v = getVoltage();
    return 0.5 * capacitance_ * v * v;
}

void Capacitor::updateState() {
    // Salva estado atual para próximo passo de tempo
    prev_voltage_ = getVoltage();
    prev_current_ = getCurrent();
}

double Capacitor::getEquivalentConductance() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Método trapezoidal: Geq = 2C/Δt
    return 2.0 * capacitance_ / timestep_;
}

double Capacitor::getEquivalentCurrent() const {
    if (timestep_ <= 0) {
        return 0.0;
    }
    // Ieq = -Geq × V(t-Δt) - I(t-Δt)
    double geq = getEquivalentConductance();
    return -geq * prev_voltage_ - prev_current_;
}

} // namespace e_XSpice