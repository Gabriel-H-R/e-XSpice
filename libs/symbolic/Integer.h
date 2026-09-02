/**
 * @file Integer.h
 * @brief Definição da classe Integer, representando números inteiros exatos.
 * 
 * Implementa números inteiros como nós terminais na árvore de expressões.
 */

#pragma once

#include "Basic.h"
#include <cstdint>

/**
 * @class Integer
 * @brief Representa um valor inteiro exato na matemática simbólica.
 */
class Integer : public Basic {
private:
    int64_t value_; ///< Valor inteiro de 64 bits armazenado.

public:
    /**
     * @brief Construtor a partir de um valor primitivo.
     * @param value O número inteiro.
     */
    explicit Integer(int64_t value);

    /** @brief Destrutor padrão */
    virtual ~Integer() = default;

    /**
     * @brief Retorna o código de tipo INTEGER.
     * @return TypeID TypeID::INTEGER.
     * @override
     */
    TypeID get_type_code() const override;

    /**
     * @brief Obtém o valor numérico armazenado.
     * @return int64_t O valor bruto.
     */
    int64_t get_value() const;

    /**
     * @brief Indica se o número é exato (Sempre verdadeiro para Integer).
     * @return bool true
     */
    bool is_exact() const;

    /**
     * @brief Indica se o número é exatamente zero.
     * @note Uma operação O(1) vital para simplificações como x + 0 = x.
     * @return bool true se for zero.
     */
    bool is_zero() const;

    /**
     * @brief Implementação do Visitor Pattern.
     * @override
     */
    void accept(Visitor &v) const override;

    /**
     * @brief Retorna os argumentos (vazio, pois é folha).
     * @override
     */
    vec_basic get_args() const override;

protected:
    std::size_t __hash__() const override;
    bool __eq__(const Basic &other) const override;
    int __cmp__(const Basic &other) const override;
};

// ============================================================================
// FUNÇÃO UTILITÁRIA PARA CRIAÇÃO RÁPIDA (Factory)
// ============================================================================
/**
 * @brief Helper para criar instâncias de Integer.
 * @param val Valor inteiro primitivo.
 * @return RCP<const Integer> O novo nó inteiro na memória.
 */
inline RCP<const Integer> integer(int64_t val) {
    return make_rcp<const Integer>(val);
}