/**
 * @file Symbol.cpp
 * @brief Implementação das lógicas de Hashing e Comparação para Symbol.
 */

#include "Symbol.h"
#include <functional> // Para std::hash

// Construtor: Inicializa a classe base implícita e move a string
Symbol::Symbol(std::string name) : name_(std::move(name)) {}

TypeID Symbol::get_type_code() const {
    return TypeID::SYMBOL;
}

const std::string& Symbol::get_name() const {
    return name_;
}

vec_basic Symbol::get_args() const {
    return {}; // Símbolos não têm filhos
}

void Symbol::accept(Visitor &v) const {
    v.visit(*this);
}

std::size_t Symbol::__hash__() const {
    // Usamos a função de hash padrão do C++ para strings.
    // É O(N) onde N é o tamanho do nome, mas lembre-se que a classe Basic
    // faz cache automático na primeira chamada, tornando chamadas futuras O(1).
    std::hash<std::string> hash_fn;
    return hash_fn(name_);
}

bool Symbol::__eq__(const Basic &other) const {
    // Fast RTTI no Basic já garantiu que 'other' é do tipo SYMBOL.
    // O static_cast aqui é 100% seguro e tem custo computacional zero!
    const Symbol &o = static_cast<const Symbol &>(other);
    return this->name_ == o.name_;
}

int Symbol::__cmp__(const Basic &other) const {
    const Symbol &o = static_cast<const Symbol &>(other);
    // string::compare retorna exatamente -1, 0 ou 1, o que é perfeito.
    int cmp = this->name_.compare(o.name_);
    if (cmp < 0) return -1;
    if (cmp > 0) return 1;
    return 0;
}