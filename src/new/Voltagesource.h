/**
 * @file VoltageSource.h
 * @brief Fonte de tensão independente
 * 
 * Uma fonte de tensão fixa a diferença de potencial entre dois nós:
 * V_plus - V_minus = V_source
 * 
 * Contribuição MNA:
 * - Adiciona uma variável extra para a corrente através da fonte
 * - Sistema aumentado: equações de tensão + equações de corrente
 * 
 * Para fonte entre nós i e j com tensão V:
 * 
 * Matriz expandida (n nós + m fontes de tensão):
 * ┌────────┬────┐ ┌───┐   ┌───┐
 * │   G    │ B  │ │ v │ = │ 0 │
 * │────────┼────│ │───│   │───│
 * │   C    │ D  │ │ i │   │ e │
 * └────────┴────┘ └───┘   └───┘
 * 
 * Onde:
 * - v: vetor de tensões nos nós
 * - i: vetor de correntes através das fontes
 * - e: vetor de tensões das fontes
 */

#ifndef VOLTAGE_SOURCE_H
#define VOLTAGE_SOURCE_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class VoltageSource : public Device {
public:
    /**
     * @brief Construtor
     * @param name Nome da fonte (ex: "V1", "VDD")
     * @param node_plus Nó positivo (potencial mais alto)
     * @param node_minus Nó negativo (potencial mais baixo)
     * @param voltage Tensão da fonte (V)
     */
    VoltageSource(const std::string& name,
                  Node* node_plus,
                  Node* node_minus,
                  double voltage);
    
    ~VoltageSource() override = default;
    
    /**
     * @brief Stamping na matriz MNA expandida
     * 
     * Fonte entre nós i e j com corrente em índice k:
     * G[i][k] += 1     (corrente sai do nó i)
     * G[j][k] -= 1     (corrente entra no nó j)
     * G[k][i] += 1     (equação de tensão)
     * G[k][j] -= 1
     * b[k] = voltage   (lado direito)
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Retorna corrente através da fonte
     * Corrente é uma incógnita extra resolvida pelo sistema
     */
    double getCurrent() const override { return current_; }
    
    /**
     * @brief Potência fornecida pela fonte P = V × I
     */
    double getPower() const override;
    
    /**
     * @brief Stamp AC: mantém a constraint KVL, mas com amplitude AC complexa.
     * A fonte de tensão AC contribui: Z[i][k]=1, Z[k][i]=1, Z[j][k]=-1, Z[k][j]=-1
     * b_ac[k] = V_ac * e^{j*phase}
     */
    void stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                 std::vector<std::complex<double>>& b_ac,
                 const std::vector<int>& node_map,
                 double omega) override;

    // Getters/Setters DC
    double getVoltage() const { return voltage_; }
    void setVoltage(double v) { voltage_ = v; }
    void setCurrent(double i) { current_ = i; }
    void setCurrentIndex(int idx) { current_index_ = idx; }
    int getCurrentIndex() const { return current_index_; }

    // Getters/Setters AC
    void setACAmplitude(double mag, double phase_deg = 0.0) {
        ac_magnitude_ = mag;
        ac_phase_rad_ = phase_deg * M_PI / 180.0;
    }
    double getACMagnitude() const { return ac_magnitude_; }
    double getACPhase() const { return ac_phase_rad_; }

private:
    double voltage_;         // Tensão DC (V)
    double current_;         // Corrente (A) — resultado da simulação
    int current_index_;      // Índice extra na matriz expandida
    double ac_magnitude_;    // Amplitude AC (V)
    double ac_phase_rad_;    // Fase AC (rad)
};

} // namespace e_XSpice

#endif // VOLTAGE_SOURCE_H