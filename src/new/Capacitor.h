/**
 * @file Capacitor.h
 * @brief Modelo de capacitor
 * 
 * Relação constitutiva: I = C × dV/dt
 * 
 * Para análise DC: Capacitor é circuito aberto (I = 0)
 * Para análise AC: Impedância Z = 1/(jωC)
 * Para análise Transiente: Usa companion model (discretização temporal)
 * 
 * Companion Model (Método de Integração Trapezoidal):
 * Em cada passo de tempo, o capacitor é substituído por:
 * - Condutância equivalente: Geq = 2C/Δt
 * - Fonte de corrente equivalente: Ieq = -Geq × V(t-Δt) - I(t-Δt)
 */

#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice {

class Capacitor : public Device {
public:
    /**
     * @brief Construtor
     * @param name Nome do capacitor
     * @param node_plus Nó positivo
     * @param node_minus Nó negativo
     * @param capacitance Capacitância (F)
     * @param initial_voltage Tensão inicial (V) para análise transiente
     */
    Capacitor(const std::string& name,
              Node* node_plus,
              Node* node_minus,
              double capacitance,
              double initial_voltage = 0.0);
    
    ~Capacitor() override = default;
    
    /**
     * @brief Stamping para análise DC (circuito aberto)
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    /**
     * @brief Stamping para análise transiente (companion model)
     * @param G Matriz
     * @param b Vetor
     * @param node_map Mapeamento de nós
     * @param timestep Passo de tempo (Δt)
     */
    void stampTransient(SparseMatrix& G, std::vector<double>& b,
                       const std::vector<int>& node_map,
                       double timestep);
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Calcula corrente através do capacitor
     * I = C × dV/dt
     */
    double getCurrent() const override;
    
    /**
     * @brief Potência instantânea P = V × I
     */
    double getPower() const override;
    
    /**
     * @brief Energia armazenada E = 0.5 × C × V²
     */
    double getEnergy() const;
    
    /**
     * @brief Atualiza estado para próximo passo de tempo
     */
    void updateState() override;
    
    /**
     * @brief Stamp AC: Y = jωC (admitância puramente imaginária)
     * Esta é a razão pela qual o capacitor bloqueia DC (ω=0 → Y=0)
     * e passa AC (ω grande → Y grande)
     */
    void stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                 std::vector<std::complex<double>>& b_ac,
                 const std::vector<int>& node_map,
                 double omega) override;

    /**
     * @brief Define passo de tempo para companion model
     */
    void setTimeStep(double dt) { timestep_ = dt; }
    
    // Getters
    double getCapacitance() const { return capacitance_; }
    double getVoltage() const;
    double getInitialVoltage() const { return initial_voltage_; }

private:
    double capacitance_;         // Capacitância (F)
    double initial_voltage_;     // Tensão inicial (V)
    double timestep_;           // Passo de tempo atual (s)
    
    // Estado anterior (para integração)
    double prev_voltage_;       // V(t-Δt)
    double prev_current_;       // I(t-Δt)
    
    /**
     * @brief Calcula condutância equivalente
     */
    double getEquivalentConductance() const;
    
    /**
     * @brief Calcula corrente equivalente (companion model)
     */
    double getEquivalentCurrent() const;
};

} // namespace e_XSpice

#endif // CAPACITOR_H