/**
 * @author Gabriel Henrique Silva
 * @file Resistor.h
 * @brief Linear resistor model
 * 
 * Ohm's law: V = I × R or I = V/R = V × G (where G = 1/R = conductance)
 * 
 * MNA contribution for resistor between nodes i and j:
 * 
 * G[i][i] += g     G[i][j] -= g
 * G[j][i] -= g     G[j][j] += g
 * 
 * where g = 1/R (conductance)
 */

#ifndef RESISTOR_H
#define RESISTOR_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class Resistor : public Device {
public:
    /**
     * @brief Constructor
     * @param name Resistor name (e.g., "R1")
     * @param node_plus Positive node
     * @param node_minus Negative node
     * @param resistance Resistance in Ohms
     */
    Resistor(const std::string& name, 
             Node* node_plus, 
             Node* node_minus,
             double resistance);
    
    /**
     * @brief Destructor
     */
    ~Resistor() override = default;
    
    /**
     * @brief Stamping in MNA matrix
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    /**
     * @brief Resistor is always linear
     */
    bool isLinear() const override { return true; }
    
    /**
     * @brief Calculates current through the resistor
     * Current flows from node_plus to node_minus
     */
    double getCurrent() const override;
    
    /**
     * @brief Calculates dissipated power P = I²R = V²/R
     */
    double getPower() const override;
    
    // Getters
    double getResistance() const { return resistance_; }
    double getConductance() const { return conductance_; }
    Node* getNodePlus() const { return nodes_[0]; }
    Node* getNodeMinus() const { return nodes_[1]; }

private:
    double resistance_;   // Resistance (Ω)
    double conductance_;  // Conductance (S) = 1/R
    
    /**
     * @brief Validates resistance
     */
    void validateResistance(double r);
};

} // namespace e_XSpice

#endif // RESISTOR_H
