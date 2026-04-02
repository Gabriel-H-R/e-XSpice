/**
 * @author Gabriel Henrique Silva
 * @file VoltageSource.h
 * @brief Independent voltage source
 * 
 * An independent voltage source fixes the potential difference between two nodes:
 * V_plus - V_minus = V_source
 * 
 * MNA contribution:
 * - Adds an extra variable for the current through the source
 * - Augmented system: voltage equations + current equations
 * 
 * For source between nodes i and j with voltage V:
 * 
 * Expanded matrix (n nodes + m voltage sources):
 * ┌────────┬────┐ ┌───┐   ┌───┐
 * │   G    │ B  │ │ v │ = │ 0 │
 * │────────┼────│ │───│   │───│
 * │   C    │ D  │ │ i │   │ e │
 * └────────┴────┘ └───┘   └───┘
 * 
 * Where:
 * - v: vector of node voltages
 * - i: vector of currents through the sources
 * - e: vector of source voltages
 */

#ifndef VOLTAGE_SOURCE_H
#define VOLTAGE_SOURCE_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class VoltageSource : public Device {
public:
    /**
     * @brief Constructor
     * @param name Source name (e.g., "V1", "VDD")
     * @param node_plus Positive node (higher potential)
     * @param node_minus Negative node (lower potential)
     * @param voltage Source voltage (V)
     */
    VoltageSource(const std::string& name,
                  Node* node_plus,
                  Node* node_minus,
                  double voltage);
    
    ~VoltageSource() override = default;
    
    /**
     * @brief Stamping in the expanded MNA matrix
     * 
     * Source between nodes i and j with current at index k:
     * G[i][k] += 1     (current leaves node i)
     * G[j][k] -= 1     (current enters node j)
     * G[k][i] += 1     (voltage equation)
     * G[k][j] -= 1
     * b[k] = voltage   (right-hand side)
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Returns current through the source
     * Current is an extra unknown solved by the system
     */
    double getCurrent() const override { return current_; }
    
    /**
     * @brief Power supplied by the source P = V × I
     */
    double getPower() const override;
    
    // Getters/Setters
    double getVoltage() const { return voltage_; }
    void setVoltage(double v) { voltage_ = v; }
    void setCurrent(double i) { current_ = i; }
    
    /**
     * @brief Sets current index in the expanded matrix
     */
    void setCurrentIndex(int idx) { current_index_ = idx; }
    int getCurrentIndex() const { return current_index_; }

private:
    double voltage_;       // Source voltage (V)
    double current_;       // Current through the source (A) - solved
    int current_index_;    // Current index in the expanded matrix
};

} // namespace e_XSpice

#endif // VOLTAGE_SOURCE_H
