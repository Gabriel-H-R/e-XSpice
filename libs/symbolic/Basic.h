/**
 * @file Basic.h
 * @brief Definição da classe Basic, a raiz do Abstract Syntax Tree (AST) matemático.
 * 
 * Todas as classes matemáticas (Symbol, Add, Mul, Integer, etc.) derivam de Basic.
 * As expressões são projetadas para serem estritamente imutáveis.
 */

#pragma once

#include "e_XSpice_rcp.h" // O nosso ponteiro super rápido
#include <vector>
#include <cstddef> // Para std::size_t

// ============================================================================
// ENUMERAÇÃO RÁPIDA DE TIPOS (Fast RTTI)
// ============================================================================
/**
 * @brief Enumeração para identificação de tipos a custo zero.
 * @note Em um CAS, checar os tipos de dados via dynamic_cast é um gargalo gigante.
 * Usamos esse TypeID para comparações ultra-rápidas durante a simplificação.
 */
enum class TypeID : unsigned int {
    SYMBOL = 0,
    INTEGER,
    RATIONAL,
    ADD,
    MUL,
    POW
    // Adicionar mais tipos conforme expandimos (Sin, Cos, Matrix, etc.)
};

// Declaração antecipada (Forward declaration) do Visitor para o Visitor Pattern
class Visitor;

// Atalho útil e onipresente no código fonte
class Basic;
using vec_basic = std::vector<RCP<const Basic>>;


// ============================================================================
// CLASSE RAIZ: BASIC
// ============================================================================
/**
 * @class Basic
 * @brief Classe base polimórfica para toda a álgebra computacional.
 * 
 * Herda de RefCounted para permitir o uso de RCP<const Basic>.
 * Aplica Hashing em Cache e ordem canônica de expressões.
 */
class Basic : public RefCounted {
private:
    /**
     * @brief Armazena o hash da expressão para evitar recalculos caros.
     * @note É `mutable` pois pode ser calculado e salvo na primeira vez 
     * que uma função constante (const) solicitar o hash. Inicialmente é 0.
     */
    mutable std::size_t hash_ = 0;

protected:
    /**
     * @brief Função pura que as classes derivadas devem implementar para gerar seu hash.
     * @return std::size_t O valor de hash específico da classe.
     */
    virtual std::size_t __hash__() const = 0;

    /**
     * @brief Compara a igualdade estrutural e matemática entre dois objetos do mesmo tipo.
     * @param other Outro objeto Basic (garantidamente do mesmo tipo pelo operador ==).
     * @return true se forem estruturalmente idênticos, false caso contrário.
     */
    virtual bool __eq__(const Basic &other) const = 0;

    /**
     * @brief Função de comparação para ordenação canônica (-1, 0, 1).
     * @note Essencial para simplificação. Ex: `y + x` é reordenado para `x + y`.
     * @param other Outro objeto Basic do mesmo tipo.
     * @return int -1 se (this < other), 0 se iguais, 1 se (this > other).
     */
    virtual int __cmp__(const Basic &other) const = 0;

public:
    /** @brief Construtor padrão e Destrutor virtual */
    Basic() = default;
    virtual ~Basic() = default;

    /**
     * @brief Retorna o código de tipo (Fast RTTI) para a classe atual.
     * @return TypeID O identificador único da classe folha.
     */
    virtual TypeID get_type_code() const = 0;

    /**
     * @brief Retorna o Hash da expressão. Usa o cache para retornar em custo O(1).
     * @return std::size_t O hash armazenado.
     */
    std::size_t hash() const noexcept;

    /**
     * @brief Compara se duas expressões são iguais (Otimizado).
     * @param other A expressão a ser comparada.
     * @return true se forem iguais.
     */
    bool eq(const Basic &other) const;

    /**
     * @brief Determina a ordem léxica/estrutural entre duas expressões.
     * @param other A expressão alvo.
     * @return int -1, 0, ou 1.
     */
    int compare(const Basic &other) const;

    /**
     * @brief Ponto de entrada para o padrão Visitor.
     * @param v Objeto visitante (Ex: impressora, derivador).
     */
    virtual void accept(Visitor &v) const = 0;

    /**
     * @brief Retorna os argumentos internos (filhos) do nó atual.
     * @note Exemplo: Para Add(x, y), retorna um vetor [x, y].
     * @return vec_basic Vetor com os filhos imutáveis (const).
     */
    virtual vec_basic get_args() const = 0;

    // --- Operadores de sobrecarga Globais ---
    friend bool operator==(const Basic &a, const Basic &b);
    friend bool operator!=(const Basic &a, const Basic &b);
};