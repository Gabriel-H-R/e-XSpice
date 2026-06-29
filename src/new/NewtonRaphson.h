/**
 * @file NewtonRaphson.h
 * @brief Solver Newton-Raphson para circuitos não-lineares
 *
 * === ALGORITMO ===
 * Para resolver F(x) = 0 (KCL com dispositivos não-lineares):
 *
 *   Inicializar x⁰ (zero ou ponto DC anterior)
 *   Para k = 0, 1, 2, ..., max_iter:
 *     Montar J(xᵏ) e F(xᵏ) via stamp() de cada dispositivo
 *     Resolver sistema linear: J·Δx = -F
 *     Atualizar: xᵏ⁺¹ = xᵏ + Δx
 *     Se ‖Δx‖ < tol: CONVERGIU
 *
 * No nosso caso, cada dispositivo chama stamp() que computa:
 *   - Parte da Jacobiana J (condutâncias linearizadas)
 *   - Parte do RHS b = -F + J*x (equivalentes Norton)
 *
 * === ESTRATÉGIAS DE CONVERGÊNCIA (inspiradas no Xyce) ===
 *
 * 1. GMIN Stepping:
 *    - Adiciona condutância GMIN de cada nó para o ground
 *    - Previne matriz singular e melhora condicionamento
 *    - Começa com GMIN grande (1e-3 S) e reduz até 1e-12 S
 *    - Cada etapa usa solução anterior como ponto de partida
 *
 * 2. Source Stepping:
 *    - Escala todas as fontes independentes por fator α ∈ [0, 1]
 *    - α=0: trivialmente V=0 converge sempre
 *    - Incrementa α em passos até α=1
 *    - Cada passo usa solução anterior como ponto de partida
 *
 * 3. PNJLIM (dentro de cada device):
 *    - Limita variação de tensão para dispositivos não-lineares
 *    - Evita overflow no cálculo de exp()
 *
 * === PARALELIZAÇÃO ===
 * - OpenMP: paraleliza montagem da Jacobiana (igual ao DC linear)
 * - MPI: distribuição de dispositivos entre processos (mesmo mecanismo do DC)
 */

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "../core/Circuit.h"
#include "../solver/SparseMatrix.h"
#include "../solver/LinearSolver.h"
#include <vector>
#include <string>
#include <functional>

namespace e_XSpice {

/**
 * @brief Status de cada tentativa de convergência
 */
enum class NRStatus {
    CONVERGED,          // Convergiu normalmente
    GMIN_CONVERGED,     // Convergiu com auxílio de GMIN stepping
    SOURCE_CONVERGED,   // Convergiu com source stepping
    FAILED              // Não convergiu com nenhuma estratégia
};

/**
 * @brief Resultado do solver Newton-Raphson
 */
struct NRResult {
    bool success;
    NRStatus status;
    int total_iterations;    // Iterações totais (inclui todos os passos)
    int nr_iterations;       // Iterações NR na tentativa final
    double residual_norm;    // Norma do resíduo final
    double elapsed_time;     // Tempo (s)
    std::string message;
    std::vector<double> solution; // Vetor solução final
};

/**
 * @brief Solver Newton-Raphson para circuitos não-lineares
 */
class NewtonRaphson {
public:
    /**
     * @brief Construtor
     * @param circuit Referência ao circuito
     */
    explicit NewtonRaphson(Circuit& circuit);

    ~NewtonRaphson() = default;

    /**
     * @brief Executa Newton-Raphson com estratégias de convergência
     *
     * Sequência de tentativas:
     *   1. NR direto (sem auxiliares)
     *   2. Se falhar: GMIN stepping
     *   3. Se falhar: Source stepping
     *
     * @param node_map    Mapeamento node_id → índice na matriz
     * @param initial_guess  Ponto de partida (vazio = usar zeros)
     * @return Resultado com status e solução
     */
    NRResult solve(const std::vector<int>& node_map,
                   const std::vector<double>& initial_guess = {});

    // Configurações
    void setMaxIterations(int n)     { max_iterations_ = n; }
    void setVTol(double tol)         { v_tol_ = tol; }   // Tolerância de tensão
    void setITol(double tol)         { i_tol_ = tol; }   // Tolerância de corrente
    void setRelTol(double tol)       { rel_tol_ = tol; }
    void setGMIN(double g)           { gmin_ = g; }
    void setVerbose(bool v)          { verbose_ = v; }
    void setLinearSolverType(SolverType t) { solver_type_ = t; }

    /**
     * @brief Define callback de progresso por iteração
     */
    void setIterationCallback(
        std::function<void(int iter, double residual)> cb) {
        iter_callback_ = cb;
    }

    /**
     * @brief Reseta o estado interno (voltages, PNJLIM)
     * Chamado antes de cada solve() independente
     */
    void reset();

private:
    Circuit& circuit_;
    int    max_iterations_;
    double v_tol_;       // Tensão: |ΔV| < v_tol (default 1µV)
    double i_tol_;       // Corrente: |ΔI| < i_tol (default 1pA)
    double rel_tol_;     // Tolerância relativa (default 1e-3)
    double gmin_;        // Condutância mínima por nó (default 1e-12 S)
    bool   verbose_;
    SolverType solver_type_;

    std::function<void(int, double)> iter_callback_;

    /**
     * @brief Uma iteração Newton-Raphson
     * Monta J(x) e b(x), resolve J·Δx = b, retorna ‖Δx‖
     */
    bool runNRIteration(const std::vector<int>& node_map,
                        double gmin_val,
                        double source_scale,
                        std::vector<double>& x,
                        double& residual);

    /**
     * @brief Monta o sistema MNA linearizado no ponto atual
     *
     * Para cada dispositivo:
     *   - Lineares: stamp normal (constante)
     *   - Não-lineares: stamp chama computação do Norton equiv. no ponto atual
     *   - GMIN: adiciona gmin_val de cada nó ao ground
     *
     * @param G          Matriz (output)
     * @param b          Vetor RHS (output)
     * @param node_map   Mapeamento
     * @param gmin_val   Condutância GMIN atual (pode ser 0)
     * @param source_scale  Fator de escala das fontes (source stepping)
     */
    void assembleSystem(SparseMatrix& G,
                        std::vector<double>& b,
                        const std::vector<int>& node_map,
                        double gmin_val,
                        double source_scale);

    /**
     * @brief Extrai solução e atualiza tensões nos nós
     */
    void extractSolution(const std::vector<double>& x,
                         const std::vector<int>& node_map);

    /**
     * @brief Verifica convergência (SPICE-style)
     * Usa ABSI + RELTOL para tensão e corrente
     */
    bool checkConvergence(const std::vector<double>& x_new,
                          const std::vector<double>& x_old,
                          int num_nodes) const;

    /**
     * @brief Tentativa de NR direto (sem auxiliares)
     */
    NRResult tryDirectNR(const std::vector<int>& node_map,
                         std::vector<double>& x);

    /**
     * @brief GMIN Stepping: reduz gmin de GMIN_START para gmin_ em passos
     */
    NRResult tryGMINStepping(const std::vector<int>& node_map,
                              std::vector<double>& x);

    /**
     * @brief Source Stepping: rampa fontes de 0 a 1 em passos
     */
    NRResult trySourceStepping(const std::vector<int>& node_map,
                                std::vector<double>& x);

    /**
     * @brief Dimensão do sistema: num_nodes + num_vsources
     */
    int systemSize(const std::vector<int>& node_map) const;
};

} // namespace e_XSpice

#endif // NEWTON_RAPHSON_H