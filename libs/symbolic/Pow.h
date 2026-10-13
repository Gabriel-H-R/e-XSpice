/**
 * @file Pow.h
 * @brief Definição da classe Pow e sua Factory canônica.
 * 
 * Representa a operação de exponenciação (Base ^ Expoente).
 * Projetada com tamanho mínimo de memória, armazenando estritamente dois ponteiros.
 */

#pragma once

#include "Basic.h"
#include "Integer.h"

// ============================================================================
// CLASSE POW (Potência)
// ============================================================================
/**
 * @class Pow
 * @brief Representação imutável de uma potência matemática.
 */
class Pow : public Basic {
private:
    RCP<const Basic> base_; ///< A base algébrica ou numérica.
    RCP<const Basic> exp_;  ///< O expoente.

public:
    /**
     * @brief Construtor estrito. Apenas a factory deve chamá-lo.
     * @param base A expressão base.
     * @param exp A expressão do expoente.
     */
    Pow(RCP<const Basic> base, RCP<const Basic> exp);

    virtual ~Pow() = default;

    TypeID get_type_code() const override;
    vec_basic get_args() const override;
    void accept(Visitor &v) const override;

    // --- Getters de Custo Zero ---
    inline const RCP<const Basic>& get_base() const noexcept { return base_; }
    inline const RCP<const Basic>& get_exp() const noexcept { return exp_; }

protected:
    std::size_t __hash__() const override;
    bool __eq__(const Basic &other) const override;
    int __cmp__(const Basic &other) const override;
};

// ============================================================================
// FACTORY CANÔNICA
// ============================================================================
/**
 * @brief Eleva uma base a um expoente aplicando regras rigorosas de simplificação.
 * @param base A expressão base.
 * @param exp O expoente.
 * @return RCP<const Basic> A expressão resultante simplificada O(1).
 */
RCP<const Basic> pow(const RCP<const Basic> &base, const RCP<const Basic> &exp);