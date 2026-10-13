/**
 * @file Mul.h
 * @brief Definição da classe Mul e sua Factory de canonicalização.
 * 
 * A classe Mul representa a multiplicação simbólica. Utiliza a mesma infraestrutura
 * de "Flattening" O(1) amortizado da classe Add, onde expressões aninhadas são 
 * esmagadas em um dicionário que mapeia as Bases para os seus Expoentes.
 */

#pragma once

#include "basic.h"
#include "integer.h"
#include "Add.h" // Reutilizamos o umap_basic_num definido no Add.h

// ============================================================================
// CLASSE MUL
// ============================================================================
/**
 * @class Mul
 * @brief Representação imutável de uma multiplicação n-ária.
 */
class Mul : public Basic {
private:
    /** @brief Coeficiente numérico multiplicativo (Ex: o '5' em '5 * x * y') */
    RCP<const Integer> coef_;
    
    /** @brief Dicionário mapeando a Base para o seu Expoente (Ex: x^2 -> Chave: x, Valor: 2) */
    umap_basic_num dict_;

public:
    /**
     * @brief Construtor estrito.
     * @param coef Coeficiente numérico constante.
     * @param dict Dicionário de bases e expoentes mesclados.
     */
    Mul(RCP<const Integer> coef, umap_basic_num dict);

    virtual ~Mul() = default;

    TypeID get_type_code() const override;
    
    /**
     * @brief Retorna os argumentos originais da multiplicação.
     * @note Como ainda não implementamos a classe Pow (Potência), 
     * retornaremos temporariamente apenas as bases. Numa versão completa, 
     * iteraríamos o dicionário retornando nós Pow(base, expoente).
     */
    vec_basic get_args() const override;
    
    void accept(Visitor &v) const override;

    // --- Getters de Custo Zero ---
    inline RCP<const Integer> get_coef() const noexcept { return coef_; }
    inline const umap_basic_num& get_dict() const noexcept { return dict_; }

protected:
    std::size_t __hash__() const override;
    bool __eq__(const Basic &other) const override;
    int __cmp__(const Basic &other) const override;
};

// ============================================================================
// FACTORY CANÔNICA
// ============================================================================
/**
 * @brief Multiplica duas expressões aplicando Flattening imediato.
 * @param a Primeiro fator.
 * @param b Segundo fator.
 * @return RCP<const Basic> A expressão resultante hiper-simplificada.
 */
RCP<const Basic> mul(const RCP<const Basic> &a, const RCP<const Basic> &b);