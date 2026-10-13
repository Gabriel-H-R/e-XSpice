/**
 * @file StrPrinter.h
 * @brief Definição do Visitor de Impressão (StrPrinter).
 * 
 * Visita os nós da Abstract Syntax Tree (AST) e constrói uma representação
 * em string legível para humanos. Usa std::ostringstream para máxima performance
 * na concatenação de caracteres.
 */

#pragma once

#include "Visitor.h"
#include "Basic.h"
#include <string>
#include <sstream>

/**
 * @class StrPrinter
 * @brief Implementação do Visitor que converte expressões matemáticas em texto.
 */
class StrPrinter : public Visitor {
private:
    /** @brief Buffer de alta performance para montagem da string. */
    std::ostringstream ss_;

public:
    StrPrinter() = default;
    virtual ~StrPrinter() = default;

    /**
     * @brief Ponto de entrada principal. Visita a árvore e retorna a string.
     * @param b Ponteiro para a raiz da expressão.
     * @return std::string A representação algébrica.
     */
    std::string apply(const RCP<const Basic> &b);

    // ========================================================================
    // MÉTODOS DE VISITA (DOUBLE DISPATCH)
    // ========================================================================
    void visit(const Symbol &x) override;
    void visit(const Integer &x) override;
    void visit(const Add &x) override;
    void visit(const Mul &x) override;
    void visit(const Pow &x) override;
};