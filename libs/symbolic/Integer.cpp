/**
 * @file Integer.cpp
 * @brief Implementação das lógicas matemáticas básicas do Integer.
 */

#include "Integer.h"
#include "Visitor.h"
#include <functional> // Para std::hash

Integer::Integer(int64_t value) : value_(value) {}

TypeID Integer::get_type_code() const {
    return TypeID::INTEGER;
}

int64_t Integer::get_value() const {
    return value_;
}

bool Integer::is_exact() const {
    return true; // Pela definição do CAS, Integers nunca são ponto-flutuante
}

bool Integer::is_zero() const {
    return value_ == 0;
}

vec_basic Integer::get_args() const {
    return {};
}

void Integer::accept(Visitor &v) const {
    v.visit(*this);
}

std::size_t Integer::__hash__() const {
    // Usamos o hash padrão do C++ para números inteiros (normalmente muito rápido)
    std::hash<int64_t> hash_fn;
    return hash_fn(value_);
}

bool Integer::__eq__(const Basic &other) const {
    // Novamente, static_cast seguro em O(1) garantido pelo TypeCode do pai
    const Integer &o = static_cast<const Integer &>(other);
    return this->value_ == o.value_;
}

int Integer::__cmp__(const Basic &other) const {
    const Integer &o = static_cast<const Integer &>(other);
    if (this->value_ < o.value_) return -1;
    if (this->value_ > o.value_) return 1;
    return 0;
}