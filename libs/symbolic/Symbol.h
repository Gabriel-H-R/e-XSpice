/**
 * @file Symbol.h
 * @brief Definição da classe Symbol, representando variáveis algébricas.
 * 
 * Símbolos são os nós folha básicos em um CAS. Identificados unicamente 
 * por seu nome (string). Como são folhas, não possuem argumentos internos.
 */

#pragma once

#include "Basic.h"
#include "Visitor.h"
#include <string>
#include <utility>

/**
 * @class Symbol
 * @brief Representa uma variável matemática contínua ou discreta (ex: 'x', 'y', 'alpha').
 */
class Symbol : public Basic {
private:
    std::string name_; ///< O nome identificador do símbolo.

public:
    /**
     * @brief Construtor principal para o Símbolo.
     * @param name A string que representa o nome da variável.
     * @note Usamos std::move para evitar cópias desnecessárias de string,
     * garantindo máxima eficiência na alocação.
     */
    explicit Symbol(std::string name);

    /**
     * @brief Destrutor virtual padrão.
     */
    virtual ~Symbol() = default;

    /**
     * @brief Retorna o código de tipo associado ao Símbolo.
     * @return TypeID O identificador rápido TypeID::SYMBOL.
     * @override
     */
    TypeID get_type_code() const override;

    /**
     * @brief Obtém o nome do símbolo.
     * @return const std::string& Referência constante para o nome.
     */
    const std::string& get_name() const;

    /**
     * @brief Implementação do Visitor Pattern.
     * @param v Objeto que visitará este nó (ex: impressora, avaliador).
     * @override
     */
    void accept(Visitor &v) const override;

    /**
     * @brief Obtém os sub-argumentos desta expressão.
     * @return vec_basic Um vetor vazio, pois Símbolos são folhas (não têm filhos).
     * @override
     */
    vec_basic get_args() const override;

protected:
    /**
     * @brief Calcula o hash estrutural do símbolo baseado em seu nome.
     * @return std::size_t Hash gerado.
     * @override
     */
    std::size_t __hash__() const override;

    /**
     * @brief Compara a igualdade entre este Símbolo e outro nó (garantido ser um Symbol).
     * @param other Outro objeto da classe Basic (que foi convertido implicitamente).
     * @return true se os nomes forem idênticos.
     * @override
     */
    bool __eq__(const Basic &other) const override;

    /**
     * @brief Compara léxicamente dois símbolos para ordenação canônica.
     * @param other Outro objeto da classe Basic.
     * @return int -1 se this < other, 0 se iguais, 1 se this > other.
     * @override
     */
    int __cmp__(const Basic &other) const override;
};

// ============================================================================
// FUNÇÃO UTILITÁRIA PARA CRIAÇÃO RÁPIDA (Factory)
// ============================================================================
/**
 * @brief Factory helper para criar símbolos de forma limpa.
 * @param name Nome do símbolo.
 * @return RCP<const Symbol> Um ponteiro inteligente apontando para o novo Símbolo.
 */
inline RCP<const Symbol> symbol(const std::string &name) {
    return make_rcp<const Symbol>(name);
}