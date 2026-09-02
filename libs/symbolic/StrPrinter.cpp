/**
 * @file StrPrinter.cpp
 * @brief Implementação das regras de formatação textual das expressões.
 */

#include "StrPrinter.h"
#include "Symbol.h"
#include "Integer.h"
#include "Add.h"
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