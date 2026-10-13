/**
 * @file Mul.cpp
 * @brief Implementação da lógica de Multiplicação e Canonicalização.
 */

#include "Mul.h"
#include "Pow.h"
#include "Visitor.h"

// ============================================================================
// MÉTODOS BÁSICOS
// ============================================================================
Mul::Mul(RCP<const Integer> coef, umap_basic_num dict)
    : coef_(std::move(coef)), dict_(std::move(dict)) {}

TypeID Mul::get_type_code() const {
    return TypeID::MUL;
}

void Mul::accept(Visitor &v) const {
    v.visit(*this);
}

vec_basic Mul::get_args() const {
    vec_basic args;
    args.reserve(dict_.size() + 1);
    
    // Na multiplicação, a identidade é 1. Só mostramos se for diferente de 1.
    if (coef_->get_value() != 1) {
        args.push_back(rcp_static_cast<const Basic>(coef_));
    }
    
    for (const auto& pair : dict_) {
        // Idealmente: args.push_back(pow(pair.first, pair.second));
        // Como não temos Pow ainda, retornamos a base.
        args.push_back(pair.first);
    }
    return args;
}

// ============================================================================
// HASH, IGUALDADE E COMPARAÇÃO
// ============================================================================
std::size_t Mul::__hash__() const {
    std::size_t h = 0;
    h ^= coef_->hash();
    for (const auto& p : dict_) {
        // A multiplicação também é comutativa, o XOR preserva o determinismo
        h ^= (p.first->hash() + p.second->hash());
    }
    return h;
}

bool Mul::__eq__(const Basic &other) const {
    const Mul &o = static_cast<const Mul &>(other);
    
    if (this->dict_.size() != o.dict_.size()) return false;
    if (!this->coef_->eq(*o.coef_)) return false;
    
    for (const auto& p : this->dict_) {
        auto it = o.dict_.find(p.first);
        if (it == o.dict_.end()) return false;
        if (!p.second->eq(*(it->second))) return false; // Verifica se o expoente é igual
    }
    return true;
}

int Mul::__cmp__(const Basic &other) const {
    const Mul &o = static_cast<const Mul &>(other);
    if (this->dict_.size() != o.dict_.size()) {
        return (this->dict_.size() < o.dict_.size()) ? -1 : 1;
    }
    return 0; 
}

// ============================================================================
// FACTORY: O MOTOR DE CANONICALIZAÇÃO (Hiper-Otimizado)
// ============================================================================
RCP<const Basic> mul(const RCP<const Basic> &a, const RCP<const Basic> &b) {
    const auto type_a = a->get_type_code();
    const auto type_b = b->get_type_code();

    // ============================================================
    // FAST PATH 1: Inteiro * Inteiro
    // ============================================================
    if (type_a == TypeID::INTEGER && type_b == TypeID::INTEGER) {
        auto ia = rcp_static_cast<const Integer>(a);
        auto ib = rcp_static_cast<const Integer>(b);
        return integer(ia->get_value() * ib->get_value());
    }

    // ============================================================
    // FAST PATH 2: Identidades do Primeiro Operando (0 e 1)
    // ============================================================
    if (type_a == TypeID::INTEGER) {
        auto ia = rcp_static_cast<const Integer>(a);
        if (ia->is_zero()) return a;            // 0 * x = 0
        if (ia->get_value() == 1) return b;     // 1 * x = x
    }

    // ============================================================
    // FAST PATH 3: Identidades do Segundo Operando (0 e 1)
    // ============================================================
    if (type_b == TypeID::INTEGER) {
        auto ib = rcp_static_cast<const Integer>(b);
        if (ib->is_zero()) return b;            // x * 0 = 0
        if (ib->get_value() == 1) return a;     // x * 1 = x
    }

    // ============================================================
    // CASO GERAL COM FLATTENING
    // ============================================================
    umap_basic_num dict;
    int64_t coef_val = 1; // Elemento neutro da multiplicação

    auto insert_term = [&](const RCP<const Basic> &term) {
        const auto type = term->get_type_code();

        // --------------------------------------------------------
        // TERMO NUMÉRICO CONSTANTE
        // --------------------------------------------------------
        if (type == TypeID::INTEGER) {
            auto integer_term = rcp_static_cast<const Integer>(term);
            coef_val *= integer_term->get_value();
            return;
        }

        // --------------------------------------------------------
        // TERMO MULTIPLICAÇÃO (O Achatamento / Flattening)
        // --------------------------------------------------------
        if (type == TypeID::MUL) {
            auto mul_node = rcp_static_cast<const Mul>(term);
            
            coef_val *= mul_node->get_coef()->get_value();

            for (const auto &[key, value] : mul_node->get_dict()) {
                auto it = dict.find(key);
                if (it != dict.end()) {
                    // Mesma base? Soma os expoentes (x^2 * x^3 = x^5)
                    const int64_t new_exp = it->second->get_value() + value->get_value();
                    
                    if (new_exp == 0) {
                        dict.erase(it); // x * x^-1 = 1 (remove do dicionário)
                    } else {
                        it->second = integer(new_exp);
                    }
                } else {
                    dict.emplace(key, value);
                }
            }
            return;
        }

        // ========================================================
        // TERMO POTÊNCIA (Pow) -> Desempacota para o dicionário!
        // ========================================================
        if (type == TypeID::POW) {
            auto pow_node = rcp_static_cast<const Pow>(term);
            
            // Verificamos se o expoente é inteiro. Nosso dicionário de Mul
            // usa Integers para manter a hiper-performance (x^2, x^3...).
            if (pow_node->get_exp()->get_type_code() == TypeID::INTEGER) {
                auto inner_base = pow_node->get_base();
                int64_t inner_exp = rcp_static_cast<const Integer>(pow_node->get_exp())->get_value();
                
                auto it = dict.find(inner_base);
                if (it != dict.end()) {
                    const int64_t new_exp = it->second->get_value() + inner_exp;
                    if (new_exp == 0) dict.erase(it);
                    else it->second = integer(new_exp);
                } else {
                    dict.emplace(inner_base, integer(inner_exp));
                }
                return; // Fim! A potência foi absorvida e achatada.
            }
        }

        // --------------------------------------------------------
        // TERMO GERAL (Símbolos, Adds, etc)
        // --------------------------------------------------------
        auto it = dict.find(term);
        if (it != dict.end()) {
            const int64_t new_exp = it->second->get_value() + 1; // Expoente base é 1
            if (new_exp == 0) {
                dict.erase(it);
            } else {
                it->second = integer(new_exp);
            }
        } else {
            dict.emplace(term, integer(1)); // Adiciona ao dicionário com expoente 1
        }
    };

    // 1. Processa o primeiro fator
    insert_term(a);

    // Short-Circuiting: Se 'a' fez o coeficiente zerar (ex: ao expandir termos),
    // cancela a avaliação da árvore 'b' e retorna instantaneamente!
    if (coef_val == 0) return integer(0);

    // 2. Processa o segundo fator
    insert_term(b);
    
    if (coef_val == 0) return integer(0);

    // ============================================================
    // CANONICALIZAÇÃO FINAL
    // ============================================================

    // 1. Todas as variáveis se anularam (Ex: x * x^-1)
    if (dict.empty()) {
        return integer(coef_val);
    }

    // 2. Sobrou apenas 1 variável com expoente 1 e coeficiente 1
    if (dict.size() == 1 && coef_val == 1) {
        auto it = dict.begin();
        if (it->second->get_value() == 1) return it->first;           // x^1 = x
        return pow(it->first, it->second);                            // Retorna o nó Pow(x, exp)!
    }

    // 3. Nó Mul genuíno
    return make_rcp<const Mul>(integer(coef_val), std::move(dict));
}