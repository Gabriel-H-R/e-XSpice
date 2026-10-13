/**
 * @file Pow.cpp
 * @brief Implementação da lógica de Exponenciação.
 */

#include "Pow.h"
#include "Visitor.h"

// ============================================================================
// MÉTODOS BÁSICOS
// ============================================================================
Pow::Pow(RCP<const Basic> base, RCP<const Basic> exp)
    : base_(std::move(base)), exp_(std::move(exp)) {}

TypeID Pow::get_type_code() const {
    return TypeID::POW;
}

void Pow::accept(Visitor &v) const {
    v.visit(*this);
}

vec_basic Pow::get_args() const {
    return {base_, exp_}; // Exatamente dois argumentos
}

// ============================================================================
// HASH, IGUALDADE E COMPARAÇÃO
// ============================================================================
std::size_t Pow::__hash__() const {
    // Como a potência NÃO é comutativa (a^b != b^a), nós não podemos usar XOR simétrico.
    // Usamos um bitshift circular seguro (<< 1) no expoente para quebrar a simetria.
    std::size_t h_base = base_->hash();
    std::size_t h_exp = exp_->hash();
    return h_base ^ (h_exp << 1);
}

bool Pow::__eq__(const Basic &other) const {
    const Pow &o = static_cast<const Pow &>(other);
    // Comparação ultrarrápida aproveitando o operator== estrutural e hashing do Basic
    return this->base_->eq(*o.base_) && this->exp_->eq(*o.exp_);
}

int Pow::__cmp__(const Basic &other) const {
    const Pow &o = static_cast<const Pow &>(other);
    // Critério: Desempata pela base. Se as bases forem iguais, desempata pelo expoente.
    int cmp_base = this->base_->compare(*o.base_);
    if (cmp_base != 0) return cmp_base;
    return this->exp_->compare(*o.exp_);
}

// ============================================================================
// FACTORY: AS LEIS DA POTÊNCIA
// ============================================================================
RCP<const Basic> pow(const RCP<const Basic> &base, const RCP<const Basic> &exp) {
    const auto type_b = base->get_type_code();
    const auto type_e = exp->get_type_code();

    // ============================================================
    // FAST PATHS COM O EXPOENTE (Identidades)
    // ============================================================
    if (type_e == TypeID::INTEGER) {
        auto ie = rcp_static_cast<const Integer>(exp);
        int64_t e_val = ie->get_value();

        if (e_val == 0) return integer(1); // x^0 = 1
        if (e_val == 1) return base;       // x^1 = x

        // ========================================================
        // COMPUTAÇÃO NUMÉRICA DIRETA (Int ^ Int)
        // ========================================================
        if (type_b == TypeID::INTEGER) {
            auto ib = rcp_static_cast<const Integer>(base);
            int64_t b_val = ib->get_value();

            if (b_val == 0 && e_val > 0) return integer(0); // 0^x = 0
            if (b_val == 1) return integer(1);              // 1^x = 1

            // Se for potência inteira positiva, calculamos na CPU instantaneamente
            if (e_val > 0) {
                int64_t result = 1;
                for (int64_t i = 0; i < e_val; ++i) {
                    result *= b_val;
                }
                return integer(result);
            }
        }

        // ========================================================
        // ACHATAMENTO DE POTÊNCIAS: (x^a)^b = x^(a * b)
        // ========================================================
        // Se a base já é um Pow e o expoente que chegou é um inteiro
        if (type_b == TypeID::POW) {
            auto inner_pow = rcp_static_cast<const Pow>(base);
            
            // Só mesclamos se o expoente interno também for Inteiro
            if (inner_pow->get_exp()->get_type_code() == TypeID::INTEGER) {
                auto inner_e = rcp_static_cast<const Integer>(inner_pow->get_exp());
                int64_t merged_exp = inner_e->get_value() * e_val;
                
                // Em vez de aninhar, retornamos uma única potência unificada
                return pow(inner_pow->get_base(), integer(merged_exp));
            }
        }
    }

    // Se nenhuma lei canônica for violada, constrói a árvore segura na memória.
    return make_rcp<const Pow>(base, exp);
}