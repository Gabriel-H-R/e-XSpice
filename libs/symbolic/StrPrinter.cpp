/**
 * @file StrPrinter.cpp
 * @brief Implementação das regras de formatação textual das expressões.
 */

#include "StrPrinter.h"
#include "Symbol.h"
#include "Integer.h"
#include "Add.h"
#include "Mul.h"
#include "Pow.h"
#include <vector>
#include <algorithm> // Para std::sort

std::string StrPrinter::apply(const RCP<const Basic> &b) {
    // Limpa o estado da stream para garantir que o Visitor 
    // possa ser reaproveitado múltiplas vezes sem recriar o objeto
    ss_.str("");
    ss_.clear();
    
    // Inicia a viagem pela árvore!
    b->accept(*this);
    
    return ss_.str();
}

void StrPrinter::visit(const Symbol &x) {
    ss_ << x.get_name();
}

void StrPrinter::visit(const Integer &x) {
    ss_ << x.get_value();
}

void StrPrinter::visit(const Add &x) {
    const auto &dict = x.get_dict();
    auto coef = x.get_coef();
    
    // 1. Extraímos os termos para podermos ordenar e garantir output determinístico
    std::vector<RCP<const Basic>> terms;
    terms.reserve(dict.size());
    for (const auto &pair : dict) {
        terms.push_back(pair.first);
    }
    
    // 2. Ordenação usando a função canônica interna de Basic (compare)
    std::sort(terms.begin(), terms.end(), [](const RCP<const Basic> &a, const RCP<const Basic> &b) {
        return a->compare(*b) == -1; // -1 significa que 'a' vem antes de 'b'
    });
    
    bool first = true;
    
    // 3. Imprime a constante primeiro (se não for zero)
    if (!coef->is_zero()) {
        coef->accept(*this);
        first = false;
    }
    
    // 4. Imprime os termos algébricos formatados
    for (const auto &base : terms) {
        if (!first) {
            ss_ << " + ";
        }
        
        auto term_coef = dict.at(base); // Busca o multiplicador no hash map
        
        // Se o coeficiente for diferente de 1, nós o imprimimos junto com o símbolo (ex: 2*x)
        if (term_coef->get_value() != 1) {
            term_coef->accept(*this);
            ss_ << "*";
        }
        
        // Imprime a base (o Símbolo)
        base->accept(*this);
        first = false;
    }
}

// ============================================================================
// IMPRESSÃO DA MULTIPLICAÇÃO (Mul)
// ============================================================================
void StrPrinter::visit(const Mul &x) {
    const auto &dict = x.get_dict();
    auto coef = x.get_coef();
    
    // 1. Extrai as bases para ordenação (garantia de determinismo)
    std::vector<RCP<const Basic>> bases;
    bases.reserve(dict.size());
    for (const auto &pair : dict) {
        bases.push_back(pair.first);
    }
    
    // 2. Ordena as bases usando a nossa função compare O(1)
    std::sort(bases.begin(), bases.end(), [](const RCP<const Basic> &a, const RCP<const Basic> &b) {
        return a->compare(*b) == -1;
    });
    
    // Flag de estado: Indica se o próximo elemento precisa de um '*' antes dele
    bool need_star = false;
    
    // 3. Imprime o coeficiente constante (se for diferente do neutro 1)
    if (coef->get_value() != 1) {
        // Formatação limpa: Se for -1, imprimimos apenas o sinal negativo se houver bases
        if (coef->get_value() == -1 && !bases.empty()) {
            ss_ << "-";
            need_star = false; // Após o '-', NÃO queremos um '*' (queremos "-x" e não "-*x")
        } else {
            coef->accept(*this);
            need_star = true;  // Se imprimiu '2', a próxima base precisa de '*' (ex: "2*x")
        }
    }
    
    // 4. Itera e imprime cada base e seu respectivo expoente
    for (const auto &base : bases) {
        
        // Adiciona o multiplicador se a flag estiver ativada
        if (need_star) {
            ss_ << "*";
        }
        
        // Verifica se a base é uma adição. Se for, precisamos de parênteses!
        bool needs_parens = (base->get_type_code() == TypeID::ADD);
        
        if (needs_parens) ss_ << "(";
        base->accept(*this); // Chama recursivamente a impressão do nó
        if (needs_parens) ss_ << ")";
        
        // 5. Verifica o expoente da base atual
        auto exp = dict.at(base);
        if (exp->get_value() != 1) {
            ss_ << "^";
            // Se o expoente for negativo, colocamos entre parênteses por elegância (ex: x^(-1))
            if (exp->get_value() < 0) {
                ss_ << "(" << exp->get_value() << ")";
            } else {
                ss_ << exp->get_value();
            }
        }
        
        // A partir do primeiro elemento algébrico, TODOS os próximos precisarão do '*'
        need_star = true;
    }
}

// ============================================================================
// POWER PRINT (Pow)
// ============================================================================
void StrPrinter::visit(const Pow &x) {
    const auto &base = x.get_base();
    const auto &exp = x.get_exp();

    // Regra Matemática Visual: Se a base for uma Adição, Multiplicação ou outra Potência,
    // precisamos cercá-la com parênteses, ex: (x + y)^2
    bool paren_base = (base->get_type_code() == TypeID::ADD || 
                       base->get_type_code() == TypeID::MUL ||
                       base->get_type_code() == TypeID::POW);

    if (paren_base) ss_ << "(";
    base->accept(*this);
    if (paren_base) ss_ << ")";

    ss_ << "^";

    // O mesmo vale para o expoente: x^(y + 1)
    bool paren_exp = (exp->get_type_code() == TypeID::ADD || 
                      exp->get_type_code() == TypeID::MUL);

    // Expoentes inteiros negativos também ganham parênteses para elegância visual: x^(-1)
    if (exp->get_type_code() == TypeID::INTEGER && rcp_static_cast<const Integer>(exp)->get_value() < 0) {
        paren_exp = true;
    }

    if (paren_exp) ss_ << "(";
    exp->accept(*this);
    if (paren_exp) ss_ << ")";
}