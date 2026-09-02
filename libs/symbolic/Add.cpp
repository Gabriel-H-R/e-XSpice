/**
 * @file Add.cpp
 * @brief Implementação da lógica de Adição e Canonicalização.
 */

#include "Add.h"
#include "Visitor.h"
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
    
    // Se houver constante diferente de zero, ela é o primeiro argumento
    if (!coef_->is_zero()) {
    // Fazemos o Upcast seguro do Integer para o Basic
    args.push_back(coef_);
}
    
    // Adiciona o restante. (Na SymEngine original as chaves são ordenadas 
    // antes de extrair os args para garantir determinismo léxico).
    for (const auto& pair : dict_) {
        // Se quiséssemos suportar multiplicação (ex: 2*x), criaríamos um 
        // nó Mul(pair.first, pair.second) aqui. Por simplicidade, retornamos a chave.
        args.push_back(pair.first);
    }
    return args;
}

// ============================================================================
// HASH, IGUALDADE E COMPARAÇÃO
// ============================================================================
std::size_t Add::__hash__() const {
    std::size_t h = 0;
    // Soma os hashes ignorando a ordem (já que a adição é comutativa).
    // O XOR (^) associado à soma embaralha os bits garantindo uma boa distribuição.
    h ^= coef_->hash();
    for (const auto& p : dict_) {
        h ^= (p.first->hash() + p.second->hash());
    }
    return h;
}

bool Add::__eq__(const Basic &other) const {
    const Add &o = static_cast<const Add &>(other);
    
    if (this->dict_.size() != o.dict_.size()) return false;
    if (!this->coef_->eq(*o.coef_)) return false;
    
    // Como é uma tabela hash, checar igualdade estrutural profunda é O(N) amortizado!
    for (const auto& p : this->dict_) {
        auto it = o.dict_.find(p.first);
        if (it == o.dict_.end()) return false;
        if (!p.second->eq(*(it->second))) return false; // Verifica o coeficiente do termo
    }
    return true;
}

int Add::__cmp__(const Basic &other) const {
    const Add &o = static_cast<const Add &>(other);
    
    // Uma heurística de ordem de expressões: primeiro pelo tamanho.
    if (this->dict_.size() != o.dict_.size()) {
        return (this->dict_.size() < o.dict_.size()) ? -1 : 1;
    }
    
    // Se necessário, faríamos uma comparação léxica ordenada termo a termo.
    // Retornamos 0 para simplificar este bloco inicial.
    return 0; 
}

// ============================================================================
// FACTORY: O MOTOR DE CANONICALIZAÇÃO
// ============================================================================
RCP<const Basic> add(const RCP<const Basic> &a, const RCP<const Basic> &b) {
    // CASO 1: Inteiro + Inteiro = Inteiro
    if (a->get_type_code() == TypeID::INTEGER && b->get_type_code() == TypeID::INTEGER) {
        auto ia = rcp_static_cast<const Integer>(a);
        auto ib = rcp_static_cast<const Integer>(b);
        return integer(ia->get_value() + ib->get_value());
    }

    // CASO 2: Expressão + 0 = Expressão
    if (a->get_type_code() == TypeID::INTEGER && rcp_static_cast<const Integer>(a)->is_zero()) return b;
    if (b->get_type_code() == TypeID::INTEGER && rcp_static_cast<const Integer>(b)->is_zero()) return a;

    // Lógica geral para 'Symbol + Symbol', 'Add + Symbol', etc.
    // Em uma implementação completa da SymEngine, nós iteraríamos sobre 'a' e 'b',
    // extraindo os inteiros para 'coef' e agrupando termos similares no 'dict'.
    
    umap_basic_num dict;
    RCP<const Integer> coef = integer(0);
    
    // Simplificação pedagógica: Assume que 'a' e 'b' são Símbolos 
    // e os armazena no dicionário com coeficiente 1.
    // (Na vida real, expandiríamos se fossem nós Add).
    
    if (a->get_type_code() == TypeID::SYMBOL) dict[a] = integer(1);
    if (b->get_type_code() == TypeID::SYMBOL) {
        if (dict.find(b) != dict.end()) {
            dict[b] = integer(dict[b]->get_value() + 1); // x + x = 2x
        } else {
            dict[b] = integer(1);
        }
    }

    return make_rcp<const Add>(coef, dict);
}