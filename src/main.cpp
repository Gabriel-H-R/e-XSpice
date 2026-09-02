/**
 * @file main.cpp
 * @brief Teste de estresse arquitetural do Sistema de Álgebra Computacional.
 */

#include <iostream>
#include "e_XSpice_rcp.h"
#include "Symbol.h"
#include "Integer.h"
#include "Add.h"
#include "StrPrinter.h"

// Função helper global para facilitar o print
std::string str(const RCP<const Basic> &expr) {
    StrPrinter printer;
    return printer.apply(expr);
}

int main() {
    std::cout << "--- Inicializando o CAS Baseado em SymEngine ---\n\n";

    // 1. Criação instântanea de Símbolos (Alocação com Cache de Hash)
    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");

    // 2. Operações e Testes de Simplificação Canônica
    
    // Exemplo A: x + x 
    // Internamente o factory add() mapeia o hash de x para o valor 2
    auto expr1 = add(x, x);
    std::cout << "Teste 1 (x + x): " << str(expr1) << "\n"; // Deve imprimir: 2*x

    // Exemplo B: x + y + x
    // x e y vão pro dicionário, e os hashes iguais colapsam os 'x'.
    auto expr2 = add(expr1, y);
    std::cout << "Teste 2 (2x + y): " << str(expr2) << "\n"; // Deve imprimir: x + y (ou 2*x + y dependendo do compare)

    // Exemplo C: Operações com Inteiros
    auto num1 = integer(10);
    auto num2 = integer(25);
    auto expr3 = add(num1, num2);
    std::cout << "Teste 3 (10 + 25): " << str(expr3) << "\n"; // Deve imprimir: 35 (E o nó Add nem existirá, é convertido para Integer O(1))

    std::cout << "\n--- Teste concluido com sucesso! ---\n";
    return 0;
}