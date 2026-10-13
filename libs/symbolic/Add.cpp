/**
 * @file Add.cpp
 * @brief Implementação da lógica de Adição e Canonicalização.
 */

#include "Add.h"
#include "Visitor.h"
#include "Mul.h"
#include <map>

// ============================================================================
// MÉTODOS BÁSICOS
// ============================================================================
Add::Add(RCP<const Integer> coef, umap_basic_num dict)
    : coef_(std::move(coef)), dict_(std::move(dict)) {}

TypeID Add::get_type_code() const {
    return TypeID::ADD;
}

void Add::accept(Visitor &v) const {
    v.visit(*this);
}

vec_basic Add::get_args() const {
    vec_basic args;
    args.reserve(dict_.size() + 1);
    
    // Se houver constante independente não-zero, ela é o primeiro argumento
    if (!coef_->is_zero()) {
        // Upcast seguro do RCP<Integer> para RCP<Basic>
        args.push_back(rcp_static_cast<const Basic>(coef_));
    }
    
    // Adiciona o restante das chaves (termos algébricos)
    for (const auto& pair : dict_) {
        args.push_back(pair.first);
    }
    return args;
}

// ============================================================================
// HASH, IGUALDADE E COMPARAÇÃO
// ============================================================================
std::size_t Add::__hash__() const {
    std::size_t h = 0;
    h ^= coef_->hash();
    for (const auto& p : dict_) {
        // O XOR simétrico das somas garante que a ordem de iteração do hash_map
        // não quebre o determinismo do hash gerado (Adição é comutativa).
        h ^= (p.first->hash() + p.second->hash());
    }
    return h;
}

bool Add::__eq__(const Basic &other) const {
    const Add &o = static_cast<const Add &>(other);
    
    // Fast path estrutural
    if (this->dict_.size() != o.dict_.size()) return false;
    if (!this->coef_->eq(*o.coef_)) return false;
    
    // Checagem profunda
    for (const auto& p : this->dict_) {
        auto it = o.dict_.find(p.first);
        if (it == o.dict_.end()) return false;
        if (!p.second->eq(*(it->second))) return false;
    }
    return true;
}

int Add::__cmp__(const Basic &other) const {
    const Add &o = static_cast<const Add &>(other);
    if (this->dict_.size() != o.dict_.size()) {
        return (this->dict_.size() < o.dict_.size()) ? -1 : 1;
    }
    // Uma implementação plena iteraria ordenadamente checando léxicamente
    return 0; 
}

// ============================================================================
// FACTORY: O MOTOR DE CANONICALIZAÇÃO (Hiper-Otimizado)
// ============================================================================
RCP<const Basic> add(const RCP<const Basic> &a, const RCP<const Basic> &b) {
    // Cache da virtual table para evitar múltiplos dispatches dinâmicos
    const auto type_a = a->get_type_code();
    const auto type_b = b->get_type_code();

    // ============================================================
    // FAST PATH 1: Inteiro + Inteiro
    // ============================================================
    if (type_a == TypeID::INTEGER && type_b == TypeID::INTEGER) {
        auto ia = rcp_static_cast<const Integer>(a);
        auto ib = rcp_static_cast<const Integer>(b);
        return integer(ia->get_value() + ib->get_value());
    }

    // ============================================================
    // FAST PATH 2: Expressão + 0 = Expressão
    // ============================================================
    if (type_a == TypeID::INTEGER) {
        auto ia = rcp_static_cast<const Integer>(a);
        if (ia->is_zero()) return b;
    }

    // ============================================================
    // FAST PATH 3: 0 + Expressão = Expressão
    // ============================================================
    if (type_b == TypeID::INTEGER) {
        auto ib = rcp_static_cast<const Integer>(b);
        if (ib->is_zero()) return a;
    }

    // ============================================================
    // CASO GERAL COM FLATTENING
    //
    // Achatamos árvores do tipo (x + y) + (z + x) diretamente 
    // para 2*x + y + z sem alocar nós Add intermediários na memória.
    // ============================================================

    // Inicializados aqui (após os fast paths) para evitar alocação de Heap
    // caso os atalhos acima encerrem a execução cedo.
    umap_basic_num dict;
    int64_t coef_val = 0;

    /**
     * @brief Processa e insere um termo algébrico ou numérico no escopo resultante.
     */
    auto insert_term = [&](const RCP<const Basic> &term) {
        const auto type = term->get_type_code();

        // --------------------------------------------------------
        // TERMO CONSTANTE NUMÉRICO
        // --------------------------------------------------------
        if (type == TypeID::INTEGER) {
            auto integer_term = rcp_static_cast<const Integer>(term);
            coef_val += integer_term->get_value();
            return;
        }

        // --------------------------------------------------------
        // TERMO ADIÇÃO (O Achatamento / Flattening)
        // --------------------------------------------------------
        if (type == TypeID::ADD) {
            auto add_node = rcp_static_cast<const Add>(term);
            
            coef_val += add_node->get_coef()->get_value();

            // Usando C++17 Structured Bindings para evitar cópias no loop
            for (const auto &[key, value] : add_node->get_dict()) {
                auto it = dict.find(key);
                if (it != dict.end()) {
                    const int64_t new_val = it->second->get_value() + value->get_value();
                    
                    // Anulação (Ex: x - x = 0)
                    if (new_val == 0) {
                        dict.erase(it);
                    } else {
                        it->second = integer(new_val);
                    }
                } else {
                    // emplace() evita inicializações redundantes do operador[]
                    dict.emplace(key, value);
                }
            }
            return;
        }

        // ========================================================
        // TERMO MULTIPLICAÇÃO
        // ========================================================
        if (type == TypeID::MUL) {
            auto mul_node = rcp_static_cast<const Mul>(term);
            const int64_t term_coef = mul_node->get_coef()->get_value();
            
            RCP<const Basic> base;
            const auto &mul_dict = mul_node->get_dict();
            
            // Otimização zero-alocação (Estado da Arte): 
            // Se o nó for apenas algo como '-1 * x' (um item no dicionário com expoente 1),
            // a base algébrica é apenas o 'x'. Nós extraímos o 'x' sem alocar nada!
            if (mul_dict.size() == 1 && mul_dict.begin()->second->get_value() == 1) {
                base = mul_dict.begin()->first;
            } else {
                // Caso contrário (ex: -5 * x^2 * y), a base é um novo nó Mul(1, {x:2, y:1})
                base = make_rcp<const Mul>(integer(1), mul_dict);
            }

            // Agora somamos o coeficiente isolado ao dicionário usando a base limpa
            auto it = dict.find(base);
            if (it != dict.end()) {
                const int64_t new_val = it->second->get_value() + term_coef;
                if (new_val == 0) {
                    dict.erase(it); // x - x = 0!
                } else {
                    it->second = integer(new_val);
                }
            } else {
                dict.emplace(base, integer(term_coef));
            }
            return;
        }

        // --------------------------------------------------------
        // TERMO GERAL (Símbolos, Funções, etc)
        // --------------------------------------------------------
        auto it = dict.find(term);
        if (it != dict.end()) {
            const int64_t new_val = it->second->get_value() + 1;
            if (new_val == 0) {
                dict.erase(it);
            } else {
                it->second = integer(new_val);
            }
        } else {
            dict.emplace(term, integer(1));
        }
    };

    // Processamento estrito e achatado de ambos os ramos
    insert_term(a);
    insert_term(b);

    // ============================================================
    // CANONICALIZAÇÃO FINAL
    // ============================================================

    // 1. Tudo se anulou -> Retorna apenas o coeficiente (Ex: 0)
    if (dict.empty()) {
        return integer(coef_val);
    }

    // 2. Sobrou apenas 1 termo algébrico (com multiplicador 1) e zero constante
    // (Ex: 0 + 1*x = x). Evita retornar Add(0, {x:1}).
    if (dict.size() == 1 && coef_val == 0) {
        auto it = dict.begin();
        if (it->second->get_value() == 1) {
            return it->first;
        }
    }

    // 3. Nó Add genuíno necessário
    return make_rcp<const Add>(integer(coef_val), std::move(dict));
}