/**
 * @file Inductor.h
 * @brief Modelo de indutor
 * 
 * Relação constitutiva: V = L × dI/dt
 * 
 * Para análise DC: Indutor é curto-circuito (V = 0)
 * Para análise AC: Impedância Z = jωL
 * Para análise Transiente: Usa companion model
 * 
 * Companion Model (Método de Integração Trapezoidal):
 * O indutor é substituído por:
 * - Resistência equivalente: Req = 2L/Δt
 * - Fonte de tensão equivalente em série
 * 
 * Como fontes de tensão adicionam variáveis, usamos Norton equivalent:
 * - Condutância: Geq = Δt/(2L)
 * - Fonte de corrente: Ieq = I(t-Δt) + Geq × V(t-Δt)
 */

#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "../core/Device.h"
#include "../core/Node.h"

namespace e_XSpice{

class Inductor : public Device {
public:
    /**
     * @brief Construtor
     * @param name Nome do indutor
     * @param node_plus Nó positivo
     * @param node_minus Nó negativo
     * @param inductance Indutância (H)
     * @param initial_current Corrente inicial (A) para análise transiente
     */
    Inductor(const std::string& name,
             Node* node_plus,
             Node* node_minus,
             double inductance,
             double initial_current = 0.0);
    
    ~Inductor() override = default;
    
    /**
     * @brief Stamping para análise DC
     * DC: indutor é curto-circuito, mas precisa de variável de corrente
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;
    
    /**
     * @brief Stamping para análise transiente (companion model)
     */
    void stampTransient(SparseMatrix& G, std::vector<double>& b,
                       const std::vector<int>& node_map,
                       double timestep);
    
    bool isLinear() const override { return true; }
    
    /**
     * @brief Retorna corrente através do indutor
     */
    double getCurrent() const override { return current_; }
    
    /**
     * @brief Potência instantânea P = V × I
     */
    double getPower() const override;
    
    /**
     * @brief Energia armazenada E = 0.5 × L × I²
     */
    double getEnergy() const;
    
    /**
     * @brief Atualiza estado para próximo passo de tempo
     */
    void updateState() override;
    
    /**
     * @brief Stamp AC: Y = 1/(jωL) = -j/(ωL)
     * Indutor bloqueia AC de alta frequência (ω grande → Y pequeno)
     * e passa DC/baixa frequência (ω→0 → curto circuito)
     * Em AC não usamos variável de corrente extra — admitância direta.
     */
    void stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                 std::vector<std::complex<double>>& b_ac,
                 const std::vector<int>& node_map,
                 double omega) override;

    /**
     * @brief Define passo de tempo
     */
    void setTimeStep(double dt) { timestep_ = dt; }
    
    /**
     * @brief Define índice da corrente na matriz (como fonte de tensão)
     */
    void setCurrentIndex(int idx) { current_index_ = idx; }
    int getCurrentIndex() const { return current_index_; }
    
    // Getters
    double getInductance() const { return inductance_; }
    double getVoltage() const;
    double getInitialCurrent() const { return initial_current_; }

private:
    double inductance_;          // Indutância (H)
    double initial_current_;     // Corrente inicial (A)
    double timestep_;           // Passo de tempo atual (s)
    double current_;            // Corrente atual (A)
    int current_index_;         // Índice na matriz expandida (para DC)
    
    // Estado anterior (para integração)
    double prev_voltage_;       // V(t-Δt)
    double prev_current_;       // I(t-Δt)
    
    /**
     * @brief Calcula condutância equivalente (Norton)
     */
    double getEquivalentConductance() const;
    
    /**
     * @brief Calcula corrente equivalente (Norton)
     */
    double getEquivalentCurrent() const;
};

} // namespace e_XSpice

#endif // INDUCTOR_H