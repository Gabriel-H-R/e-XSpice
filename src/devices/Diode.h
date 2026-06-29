/**
 * @file Diode.h
 * @brief Modelo de diodo com equação de Shockley e PNJLIM
 *
 * Inspirado no modelo D do Xyce/SPICE.
 *
 * === MODELO FÍSICO ===
 * Equação de Shockley (diodo ideal):
 *
 *   I_D = I_S * (exp(V_D / (n * V_T)) - 1)
 *
 * Onde:
 *   I_S  = corrente de saturação reversa  (típico: 1e-14 A)
 *   n    = fator de idealidade            (1 ≤ n ≤ 2)
 *   V_T  = tensão térmica = kT/q          (≈26mV a 300K)
 *   V_D  = tensão anode-catodo
 *
 * === LINEARIZAÇÃO PARA NEWTON-RAPHSON ===
 * A cada iteração k, calcula o equivalente Norton:
 *
 *   G_eq = dI_D/dV|_{V_k} = (I_S / (n*V_T)) * exp(V_k/(n*V_T))
 *                         = (I_D + I_S) / (n*V_T)   [mais estável]
 *   I_eq = I_D(V_k) - G_eq * V_k   [fonte de corrente equivalente]
 *
 * Stamp: adiciona G_eq como condutância e I_eq como fonte de corrente.
 * Isso é idêntico ao que Xyce faz internamente.
 *
 * === PNJLIM — VOLTAGE LIMITING ===
 * Sem limitação, exp(V/nVt) causa overflow numérico em poucas iterações.
 * O PNJLIM (PN Junction LIMiter) restringe a variação de tensão por passo:
 *
 *   Se V_new > V_old + n*V_T:
 *     V_lim = V_old + n*V_T * ln((V_new - V_old)/(n*V_T) + 1)
 *   Se V_new < -5*n*V_T:
 *     V_lim = -5*n*V_T
 *
 * === ANÁLISE AC (SMALL-SIGNAL) ===
 * No ponto DC V_Q:
 *   G_d = dI_D/dV|_{V_Q} = (I_D + I_S) / (n*V_T)
 *
 * O diodo vira uma condutância G_d em paralelo com a junção (real pura).
 * A capacitância de junção C_j é a contribuição imaginária em AC.
 */

#ifndef DIODE_H
#define DIODE_H

#include "../core/Device.h"
#include "../core/Node.h"
#include <string>
#include <vector>
#include <complex>
#include <cmath>

namespace e_XSpice {


struct DiodeModelParams {
    double IS  = 1e-14;   // Corrente de saturação (A)
    double N   = 1.0;     // Fator de idealidade
    double RS  = 0.0;     // Resistência série (Ω)
    double CJO = 0.0;     // Capacitância de junção zero-bias (F)
    double VJ  = 0.7;     // Potencial de junção (V)
    double M   = 0.5;     // Gradiente de junção
    double T   = 300.15;  // Temperatura (K)
};

class Diode : public Device {
public:
    using ModelParams = DiodeModelParams;
    
    /**
     * @brief Construtor
     * @param name    Nome do diodo (ex: "D1")
     * @param anode   Nó anodo (positivo)
     * @param cathode Nó catodo (negativo)
     * @param params  Parâmetros do modelo (opcional, usa padrões SPICE)
     */
    Diode(const std::string& name,
          Node* anode,
          Node* cathode,
          const ModelParams& params = ModelParams());

    ~Diode() override = default;

    /**
     * @brief Stamping Newton-Raphson: equivalente Norton linearizado
     *
     * Lê V_D atual dos nós, aplica PNJLIM, calcula G_eq e I_eq.
     * Chamado em cada iteração NR — por isso o stamp é não-linear.
     */
    void stamp(SparseMatrix& G, std::vector<double>& b,
               const std::vector<int>& node_map) override;

    /**
     * @brief Stamping AC: small-signal conductance G_d no ponto DC
     * Opcionalmente inclui capacitância de junção C_j.
     */
    void stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                 std::vector<std::complex<double>>& b_ac,
                 const std::vector<int>& node_map,
                 double omega) override;

    bool isLinear() const override { return false; }

    /**
     * @brief Corrente DC através do diodo (Shockley)
     */
    double getCurrent() const override;

    /**
     * @brief Potência instantânea
     */
    double getPower() const override;

    /**
     * @brief Tensão ânodo-cátodo
     */
    double getVoltage() const;

    /**
     * @brief Condutância de small-signal no ponto operação atual
     */
    double getSmallSignalConductance() const;

    /**
     * @brief Capacitância de junção (para AC)
     */
    double getJunctionCapacitance() const;

    const ModelParams& getParams() const { return params_; }

private:
    ModelParams params_;
    double vd_prev_;     // Tensão V_D da iteração anterior (para PNJLIM)
    bool   first_iter_;  // Flag para primeira iteração

    /**
     * @brief Tensão térmica V_T = kT/q
     */
    double thermalVoltage() const;

    /**
     * @brief Corrente Shockley para tensão V
     */
    double shockleyCurrent(double V) const;

    /**
     * @brief Condutância (derivada da corrente) para tensão V
     */
    double shockleyConductance(double V) const;

    /**
     * @brief PNJLIM — limita variação de V para evitar overflow exp()
     * Algoritmo direto do SPICE3 (pnjlim.c)
     *
     * @param v_new  Nova tensão (solução atual)
     * @param v_old  Tensão anterior (iteração k-1)
     * @param vt     Tensão térmica n*V_T
     * @param vcrit  Tensão crítica (threshold para limitar)
     * @return Tensão limitada
     */
    double pnjlim(double v_new, double v_old, double vt, double vcrit) const;
};

} // namespace e_XSpice

#endif // DIODE_H