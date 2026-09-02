/**
 * @file Basic.c
 * @brief Implementação dos métodos não-virtuais da classe Basic.
 */

#include "Basic.h"
#include <typeinfo> // Usado apenas fallback defensivo, evitado para performance

// ============================================================================
// HASHING E CACHE
// ============================================================================
std::size_t Basic::hash() const noexcept {
    // Checa se o hash já foi computado e salvo em cache (diferente de zero).
    if (hash_ == 0) {
        // Se for 0, chama o gerador de hash específico da classe derivada.
        hash_ = this->__hash__();
        
        // Garantia de anti-colisão com o valor zero:
        // Se a classe gerar 0 por azar, forçamos a ser 1 para evitar recalcular 
        // a mesma árvore no futuro.
        if (hash_ == 0) {
            hash_ = 1;
        }
    }
    return hash_;
}

// ============================================================================
// COMPARAÇÃO DE IGUALDADE (eq)
// ============================================================================
bool Basic::eq(const Basic &other) const {
    // Verificação de mesma identidade de memória (O ponteiro é o mesmo?).
    // Se eles apontam pro mesmo endereço, são obviamente iguais (O(1)).
    if (this == &other) {
        return true;
    }

    // O Hash é diferente?
    // Se o hash for diferente, as expressões obrigatoriamente são diferentes.
    // Isso evita navegar numa árvore gigantesca desnecessariamente.
    if (this->hash() != other.hash()) {
        return false;
    }

    // Tipos diferentes?
    // Se um é Symbol e o outro é Add, não podem ser iguais.
    if (this->get_type_code() != other.get_type_code()) {
        return false;
    }

    // Se passou por todas as barreiras rápidas acima (mesmo hash, mesmo tipo,
    // mas memórias diferentes), caímos na checagem profunda e segura específica da classe.
    return this->__eq__(other);
}

// ============================================================================
// ORDENAÇÃO CANÔNICA (compare)
// ============================================================================
int Basic::compare(const Basic &other) const {
    // Se forem a mesma instância na memória
    if (this == &other) return 0;

    // Fast RTTI: Os tipos são diferentes?
    TypeID id1 = this->get_type_code();
    TypeID id2 = other.get_type_code();

    // Uma regra rigorosa da álgebra computacional é a ordenação global de tipos.
    // Exemplo: Symbols sempre vêm antes de Funções.
    if (id1 != id2) {
        return (id1 < id2) ? -1 : 1;
    }

    // Se os tipos forem iguais e os endereços diferentes, o desempate
    // profundo é feito pela implementação específica (ex: 'x' comparado com 'y').
    return this->__cmp__(other);
}

// ============================================================================
// SOBRECARGA DOS OPERADORES == E !=
// ============================================================================
bool operator==(const Basic &a, const Basic &b) {
    return a.eq(b);
}

bool operator!=(const Basic &a, const Basic &b) {
    return !a.eq(b);
}