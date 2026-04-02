/**
 * @author Gabriel Henrique Silva
 * @file CurrentSource.h
 * @brief Independent current source
 * 
 * A current source injects a fixed current into the circuit.
 * Convention: current flows from the negative node to the positive.
 * 
 * MNA contribution for current source I between nodes i and j:
 * - Does not add extra rows/columns (unlike voltage source)
 * - Only modifies vector b:
 * 
 * b[i] -= I  (current leaves node i)
 * b[j] += I  (current enters node j)
 * 
 * Example: I = 1A from node 0 (GND) to node 1
 * → b[1] += 1.0 (1A enters node 1)
 */

#ifndef CURRENT_SOURCE_H
#define CURRENT_SOURCE_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class CurrentSource : public Device {
public:
    /**
     * @brief Constructor
     * @param name Source name (e.g., "I1", "Ibias")
     * @param node_plus Node where current flows into (input)
     * @param node_minus Node where current flows out (output)
     * @param current Source current (A)
     */
    CurrentSource(const std::string& name,
                  Node* node_plus,
                  Node* node_minus,
                  double current);
    
    ~CurrentSource() override = default;
    
    /**
     * @brief Stamping in MNA matrix
     * Current source only modifies vector b
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Returns source current (known)
     */
    double getCurrent() const override { return current_; }
    
    /**
     * @brief Power supplied P = V × I
     * where V = V_plus - V_minus
     */
    double getPower() const override;
    
    // Getters/Setters
    double getCurrentValue() const { return current_; }
    void setCurrent(double i) { current_ = i; }

private:
    double current_;  // Source current (A)
};

} // namespace e_XSpice

#endif // CURRENT_SOURCE_H
