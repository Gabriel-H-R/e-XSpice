/**
 * @file Resistor.h
 * @brief Modelo de resistor linear
 * 
 * Lei de Ohm: V = I × R ou I = V/R = V × G (onde G = 1/R = condutância)
 * 
 * Contribuição MNA para resistor entre nós i e j:
 * 
 * G[i][i] += g     G[i][j] -= g
 * G[j][i] -= g     G[j][j] += g
 * 
 * onde g = 1/R (condutância)
 */

#ifndef RESISTOR_H
#define RESISTOR_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class Resistor : public Device {
public:
    /**
     * @brief Construtor
     * @param name Nome do resistor (ex: "R1")
     * @param node_plus Nó positivo
     * @param node_minus Nó negativo
     * @param resistance Resistência em Ohms
     */
    Resistor(const std::string& name, 
             Node* node_plus, 
             Node* node_minus,
             double resistance);
    
    /**
     * @brief Destrutor
     */
    ~Resistor() override = default;
    
    /**
     * @brief Stamping na matriz MNA
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    /**
     * @brief Resistor é sempre linear
     */
    bool isLinear() const override { return true; }
    
    /**
     * @brief Calcula corrente através do resistor
     * Corrente vai de node_plus para node_minus
     */
    double getCurrent() const override;
    
    /**
     * @brief Calcula potência dissipada P = I²R = V²/R
     */
    double getPower() const override;
    
    /**
     * @brief Stamp AC: resistor contribui Y = 1/R (real, sem variação com ω)
     */
    void stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                 std::vector<std::complex<double>>& b_ac,
                 const std::vector<int>& node_map,
                 double omega) override;

    // Getters
    double getResistance() const { return resistance_; }
    double getConductance() const { return conductance_; }
    Node* getNodePlus() const { return nodes_[0]; }
    Node* getNodeMinus() const { return nodes_[1]; }

private:
    double resistance_;   // Resistência (Ω)
    double conductance_;  // Condutância (S) = 1/R
    
    /**
     * @brief Valida resistência
     */
    void validateResistance(double r);
};

} // namespace e_XSpice

#endif // RESISTOR_H