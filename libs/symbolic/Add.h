/**
 * @file Add.h
 * @brief Definição da classe Add, representando operações de soma (+).
 * 
 * A classe Add achata as somas e usa Dicionários Hash para combinar 
 * termos semelhantes instantaneamente em tempo O(1).
 */

#pragma once

#include "Basic.h"
#include "Integer.h"
#include <unordered_map>
#include <map>

// ============================================================================
// HASHING CUSTOMIZADO PARA MAPAS (O Segredo da Performance)
// ============================================================================
/**
 * @brief Functor para extrair o hash de um RCP<const Basic>.
 * @note Delegamos diretamente para o método O(1) hash() da classe Basic.
 */
struct RCPBasicHash {
    std::size_t operator()(const RCP<const Basic>& a) const noexcept {
        return a->hash();
    }
};

/**
 * @brief Functor para checar igualdade de chaves no mapa.
 * @note Delegamos para o método ultra-rápido eq() de Basic.
 */
struct RCPBasicEq {
    bool operator()(const RCP<const Basic>& a, const RCP<const Basic>& b) const {
        return a->eq(*b);
    }
};

/** 
 * @brief Dicionário fundamental para Otimização de Comutatividade.
 * Mapeia a Base (Símbolo/Expressão) para o seu Coeficiente Numérico.
 * Ex: 3x -> Chave: 'x', Valor: Integer(3).
 */
using umap_basic_num = std::unordered_map<RCP<const Basic>, RCP<const Integer>, RCPBasicHash, RCPBasicEq>;


// ============================================================================
// CLASSE ADD
// ============================================================================
/**
 * @class Add
 * @brief Representa uma soma de N termos matemáticos.
 */
class Add : public Basic {
private:
    /** @brief Termo constante numérico (Ex: o '5' em 'x + y + 5'). */
    RCP<const Integer> coef_;
    
    /** @brief Dicionário de termos e seus coeficientes. Ex: 'x' -> 2, 'y' -> 1. */
    umap_basic_num dict_;

public:
    /**
     * @brief Construtor estrito.
     * @note De acordo com o design da SymEngine, construtores são geralmente 
     * ocultos/protegidos, mas deixamos público para facilitar a factory. 
     * Nunca instancie via 'new Add(...)', use sempre a função 'add()'.
     * @param coef Coeficiente numérico exato.
     * @param dict Dicionário de termos algébricos.
     */
    Add(RCP<const Integer> coef, umap_basic_num dict);

    virtual ~Add() = default;

    TypeID get_type_code() const override;
    vec_basic get_args() const override;
    void accept(Visitor &v) const override;

    /**
     * @brief Obtém o coeficiente numérico constante da adição.
     * @return RCP<const Integer> O termo constante.
     */
    inline RCP<const Integer> get_coef() const noexcept { return coef_; }

    /**
     * @brief Obtém o dicionário interno de termos e coeficientes.
     * @return const umap_basic_num& Referência const para o mapa hash.
     */
    inline const umap_basic_num& get_dict() const noexcept { return dict_; }

protected:
    std::size_t __hash__() const override;
    bool __eq__(const Basic &other) const override;
    int __cmp__(const Basic &other) const override;
};

// ============================================================================
// FUNÇÃO FACTORY (Construtor Canônico)
// ============================================================================
/**
 * @brief Cria uma adição algébrica simplificando-a o máximo possível.
 * @param a Primeiro termo.
 * @param b Segundo termo.
 * @return RCP<const Basic> A expressão canônica resultante (que pode ser rebaixada a Integer ou Symbol).
 */
RCP<const Basic> add(const RCP<const Basic> &a, const RCP<const Basic> &b);