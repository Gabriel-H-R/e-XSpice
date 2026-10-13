/**
 * @file Visitor.h
 * @brief Definição da interface Visitor para percorrimento do AST (Abstract Syntax Tree).
 * 
 * Seguindo o design da SymEngine, algoritmos como Derivação, Substituição e Impressão 
 * não são implementados dentro das classes matemáticas. Eles são implementados como 
 * Visitors que viajam pelos nós da árvore.
 */

#pragma once

// Forward declarations para evitar dependências circulares pesadas
class Symbol;
class Integer;
class Add;
class Mul;
class Pow;

/**
 * @class Visitor
 * @brief Interface base para o padrão Visitor.
 * 
 * Cada classe concreta derivada de Basic invocará o método visit() correspondente
 * ao seu próprio tipo, garantindo o "Double Dispatch" e o RTTI seguro.
 */
class Visitor {
public:
    virtual ~Visitor() = default;

    /**
     * @brief Visita um nó do tipo Símbolo.
     * @param x O Símbolo a ser processado.
     */
    virtual void visit(const Symbol &x) = 0;

    /**
     * @brief Visita um nó do tipo Inteiro.
     * @param x O Inteiro a ser processado.
     */
    virtual void visit(const Integer &x) = 0;

    /**
     * @brief Visita um nó do tipo Adição.
     * @param x A Adição a ser processada.
     */
    virtual void visit(const Add &x) = 0;

    /**
     * @brief Visita um nó do tipo Multiplicação.
     * @param x A Multiplicação a ser processada.
     */
    virtual void visit(const Mul &x) = 0;

    /**
     * @brief Visita um nó do tipo Potência.
     * @param x A Potência a ser processada.
     */
    virtual void visit(const Pow &x) = 0;
};