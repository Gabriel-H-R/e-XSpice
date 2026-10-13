/**
 * @file main.cpp
 * @brief Suíte de testes do CAS.
 *
 * Testa as funcionalidades atualmente disponíveis:
 *   - Symbol
 *   - Integer
 *   - Add
 *   - Mul
 *   - Subtração via Add + Mul
 *   - Canonicalização
 *   - Combinação de termos
 *   - Combinação de fatores
 *   - Identidade 0 / 1
 *   - Comutatividade
 *   - Associatividade
 *   - Cancelamento
 *   - Expressões aninhadas
 */

#include <iostream>
#include <string>
#include <stdio.h>
#include <memory>
#include <fstream>

#include "../libs/symbolic/e_XSpice_rcp.h"
#include "../libs/symbolic/Symbol.h"
#include "../libs/symbolic/Integer.h"
#include "../libs/symbolic/Add.h"
#include "../libs/symbolic/Mul.h"
#include "../libs/symbolic/Pow.h"
#include "../libs/symbolic/StrPrinter.h"


// ============================================================
// Helpers
// ============================================================

std::string str(const RCP<const Basic>& expr)
{
    StrPrinter printer;
    return printer.apply(expr);
}


RCP<const Basic> sub(
    const RCP<const Basic>& a,
    const RCP<const Basic>& b)
{
    // a - b = a + (-1 * b)

    auto minus_one = integer(-1);

    return add(a, mul(minus_one, b));
}

void expect_equal(
    const RCP<const Basic>& a,
    const RCP<const Basic>& b,
    const std::string& test_name
) {
    if (*a == *b) {
        std::cout << "[PASS] " << test_name << "\n";
    } else {
        std::cout << "[FAIL] " << test_name << "\n";
        std::cout << "       A = " << str(a) << "\n";
        std::cout << "       B = " << str(b) << "\n";
    }
}


// ============================================================
// Sistema simples de testes
// ============================================================

int tests_run = 0;
int tests_passed = 0;
int tests_failed = 0;


void test(
    const std::string& name,
    const RCP<const Basic>& expression,
    const std::string& expected)
{
    tests_run++;

    std::string result = str(expression);

    if (result == expected)
    {
        tests_passed++;

        std::cout
            << "[PASS] "
            << name
            << "\n"
            << "       Resultado: " << result
            << "\n";
    }
    else
    {
        tests_failed++;

        std::cout
            << "[FAIL] "
            << name
            << "\n"
            << "       Esperado : " << expected
            << "\n"
            << "       Obtido   : " << result
            << "\n";
    }

    std::cout << "\n";
}


// ============================================================
// MAIN
// ============================================================

int main()
{
    std::cout
        << "============================================\n"
        << "      e-XSpice - CAS Test Suite\n"
        << "============================================\n\n";


    // ========================================================
    // Símbolos básicos
    // ========================================================

    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");


    std::cout
        << "============================================\n"
        << "1. SYMBOLS\n"
        << "============================================\n\n";


    test(
        "Símbolo x",
        x,
        "x"
    );

    test(
        "Símbolo y",
        y,
        "y"
    );


    // ========================================================
    // Inteiros
    // ========================================================

    std::cout
        << "============================================\n"
        << "2. INTEGERS\n"
        << "============================================\n\n";


    test(
        "Integer 0",
        integer(0),
        "0"
    );

    test(
        "Integer 1",
        integer(1),
        "1"
    );

    test(
        "Integer 2",
        integer(2),
        "2"
    );

    test(
        "Integer -1",
        integer(-1),
        "-1"
    );

    test(
        "Integer 42",
        integer(42),
        "42"
    );


    // ========================================================
    // ADD básico
    // ========================================================

    std::cout
        << "============================================\n"
        << "3. ADD - BASIC\n"
        << "============================================\n\n";


    test(
        "x + y",
        add(x, y),
        "x + y"
    );

    test(
        "x + x",
        add(x, x),
        "2*x"
    );

    test(
        "y + y",
        add(y, y),
        "2*y"
    );


    // ========================================================
    // Identidade aditiva
    // ========================================================

    std::cout
        << "============================================\n"
        << "4. ADD - IDENTITY\n"
        << "============================================\n\n";


    test(
        "x + 0",
        add(x, integer(0)),
        "x"
    );

    test(
        "0 + x",
        add(integer(0), x),
        "x"
    );

    test(
        "0 + 0",
        add(integer(0), integer(0)),
        "0"
    );

    test(
        "5 + 0",
        add(integer(5), integer(0)),
        "5"
    );


    // ========================================================
    // Soma de inteiros
    // ========================================================

    std::cout
        << "============================================\n"
        << "5. ADD - INTEGERS\n"
        << "============================================\n\n";


    test(
        "2 + 3",
        add(integer(2), integer(3)),
        "5"
    );

    test(
        "10 + 25",
        add(integer(10), integer(25)),
        "35"
    );

    test(
        "10 + (-10)",
        add(integer(10), integer(-10)),
        "0"
    );


    // ========================================================
    // Comutatividade do ADD
    // ========================================================

    std::cout
        << "============================================\n"
        << "6. ADD - COMMUTATIVITY\n"
        << "============================================\n\n";


    test(
        "x + y",
        add(x, y),
        str(add(y, x))
    );

    test(
        "x + 2*x",
        add(x, mul(integer(2), x)),
        str(add(mul(integer(2), x), x))
    );


    // ========================================================
    // Associatividade do ADD
    // ========================================================

    std::cout
        << "============================================\n"
        << "7. ADD - ASSOCIATIVITY\n"
        << "============================================\n\n";


    auto add_left =
        add(
            add(x, y),
            z
        );

    auto add_right =
        add(
            x,
            add(y, z)
        );


    test(
        "(x + y) + z",
        add_left,
        str(add_right)
    );


    // ========================================================
    // Combinação de termos
    // ========================================================

    std::cout
        << "============================================\n"
        << "8. ADD - LIKE TERMS\n"
        << "============================================\n\n";


    test(
        "x + x + x",
        add(add(x, x), x),
        "3*x"
    );

    test(
        "2*x + x",
        add(mul(integer(2), x), x),
        "3*x"
    );

    test(
        "x + 2*x",
        add(x, mul(integer(2), x)),
        "3*x"
    );

    test(
        "2*x + 3*x",
        add(
            mul(integer(2), x),
            mul(integer(3), x)
        ),
        "5*x"
    );


    // ========================================================
    // SUBTRAÇÃO
    // ========================================================

    std::cout
        << "============================================\n"
        << "9. SUBTRACTION\n"
        << "============================================\n\n";


    test(
        "x - x",
        sub(x, x),
        "0"
    );

    test(
        "x - 0",
        sub(x, integer(0)),
        "x"
    );

    test(
        "0 - x",
        sub(integer(0), x),
        "-x"
    );

    test(
        "5 - 5",
        sub(integer(5), integer(5)),
        "0"
    );

    test(
        "10 - 3",
        sub(integer(10), integer(3)),
        "7"
    );


    // ========================================================
    // MUL básico
    // ========================================================

    std::cout
        << "============================================\n"
        << "10. MUL - BASIC\n"
        << "============================================\n\n";


    test(
        "x * y",
        mul(x, y),
        "x*y"
    );

    test(
        "x * x",
        mul(x, x),
        "x^2"
    );

    test(
        "y * y",
        mul(y, y),
        "y^2"
    );


    // ========================================================
    // Identidade multiplicativa
    // ========================================================

    std::cout
        << "============================================\n"
        << "11. MUL - IDENTITY\n"
        << "============================================\n\n";


    test(
        "x * 1",
        mul(x, integer(1)),
        "x"
    );

    test(
        "1 * x",
        mul(integer(1), x),
        "x"
    );

    test(
        "x * 0",
        mul(x, integer(0)),
        "0"
    );

    test(
        "0 * x",
        mul(integer(0), x),
        "0"
    );

    test(
        "0 * 0",
        mul(integer(0), integer(0)),
        "0"
    );


    // ========================================================
    // Multiplicação por inteiros
    // ========================================================

    std::cout
        << "============================================\n"
        << "12. MUL - INTEGERS\n"
        << "============================================\n\n";


    test(
        "2 * 3",
        mul(integer(2), integer(3)),
        "6"
    );

    test(
        "10 * 5",
        mul(integer(10), integer(5)),
        "50"
    );

    test(
        "-1 * x",
        mul(integer(-1), x),
        "-x"
    );

    test(
        "2 * x",
        mul(integer(2), x),
        "2*x"
    );


    // ========================================================
    // Comutatividade do MUL
    // ========================================================

    std::cout
        << "============================================\n"
        << "13. MUL - COMMUTATIVITY\n"
        << "============================================\n\n";


    test(
        "x * y",
        mul(x, y),
        str(mul(y, x))
    );

    test(
        "2 * x",
        mul(integer(2), x),
        str(mul(x, integer(2)))
    );


    // ========================================================
    // Associatividade do MUL
    // ========================================================

    std::cout
        << "============================================\n"
        << "14. MUL - ASSOCIATIVITY\n"
        << "============================================\n\n";


    auto mul_left =
        mul(
            mul(x, y),
            z
        );

    auto mul_right =
        mul(
            x,
            mul(y, z)
        );


    test(
        "(x * y) * z",
        mul_left,
        str(mul_right)
    );


    // ========================================================
    // Combinação de fatores
    // ========================================================

    std::cout
        << "============================================\n"
        << "15. MUL - FACTOR COMBINATION\n"
        << "============================================\n\n";


    test(
        "x^2 * x",
        mul(
            mul(x, x),
            x
        ),
        "x^3"
    );

    test(
        "x^2 * x^2",
        mul(
            mul(x, x),
            mul(x, x)
        ),
        "x^4"
    );

    test(
        "x^2 * x * y",
        mul(
            mul(x, x),
            mul(x, y)
        ),
        "x^3*y"
    );


    // ========================================================
    // Múltiplos símbolos
    // ========================================================

    std::cout
        << "============================================\n"
        << "16. MULTIPLE SYMBOLS\n"
        << "============================================\n\n";


    test(
        "x * y * z",
        mul(
            mul(x, y),
            z
        ),
        "x*y*z"
    );

    test(
        "x * x * y * y",
        mul(
            mul(x, x),
            mul(y, y)
        ),
        "x^2*y^2"
    );

    test(
        "x^2 * y^2 * x",
        mul(
            mul(
                mul(x, x),
                mul(y, y)
            ),
            x
        ),
        "x^3*y^2"
    );


    // ========================================================
    // Expressões aninhadas
    // ========================================================

    std::cout
        << "============================================\n"
        << "17. NESTED EXPRESSIONS\n"
        << "============================================\n\n";


    test(
        "2 * (x + y)",
        mul(
            integer(2),
            add(x, y)
        ),
        "2*(x + y)"
    );

    test(
        "(x + y) * 2",
        mul(
            add(x, y),
            integer(2)
        ),
        "2*(x + y)"
    );


    // ========================================================
    // Soma dentro de soma
    // ========================================================

    std::cout
        << "============================================\n"
        << "18. NESTED ADD\n"
        << "============================================\n\n";


    test(
        "(x + y) + z",
        add(
            add(x, y),
            z
        ),
        str(
            add(
                x,
                add(y, z)
            )
        )
    );


    // ========================================================
    // Expressões mais complexas
    // ========================================================

    std::cout
        << "============================================\n"
        << "19. COMPLEX EXPRESSIONS\n"
        << "============================================\n\n";


    auto expr1 =
        mul(
            add(x, y),
            integer(2)
        );


    test(
        "2 * (x + y)",
        expr1,
        "2*(x + y)"
    );


    auto expr2 =
        add(
            mul(integer(2), x),
            mul(integer(3), x)
        );


    test(
        "2*x + 3*x",
        expr2,
        "5*x"
    );


    auto expr3 =
        sub(
            mul(integer(3), x),
            mul(integer(3), x)
        );


    test(
        "3*x - 3*x",
        expr3,
        "0"
    );


    auto expr4 =
        mul(
            add(x, y),
            add(x, y)
        );


    test(
        "(x + y) * (x + y)",
        expr4,
        "(x + y)^2"
    );


    // ========================================================
    // Testes de combinação mais agressivos
    // ========================================================

    std::cout
        << "============================================\n"
        << "20. STRESS - CANONICALIZATION\n"
        << "============================================\n\n";


    auto stress1 =
        mul(
            x,
            mul(
                y,
                mul(x, y)
            )
        );


    test(
        "x * y * x * y",
        stress1,
        "x^2*y^2"
    );


    auto stress2 =
        add(
            add(x, y),
            add(x, y)
        );


    test(
        "(x + y) + (x + y)",
        stress2,
        "2*x + 2*y"
    );


    auto stress3 =
        add(
            add(
                mul(integer(2), x),
                mul(integer(3), y)
            ),
            add(
                mul(integer(4), x),
                mul(integer(5), y)
            )
        );


    test(
        "2*x + 3*y + 4*x + 5*y",
        stress3,
        "6*x + 8*y"
    );

    auto stress4 =
        mul(
            integer(2),
            add(x, y)
        );


    test(
        "2*(x + y)",
        stress4,
        "2*(x + y)"
    );


    // ========================================================
    // Testes de cancelamento
    // ========================================================

    std::cout
        << "============================================\n"
        << "21. CANCELLATION\n"
        << "============================================\n\n";


    test(
        "x - x",
        sub(x, x),
        "0"
    );

    test(
        "2*x - 2*x",
        sub(
            mul(integer(2), x),
            mul(integer(2), x)
        ),
        "0"
    );

    test(
        "3*x + (-3*x)",
        add(
            mul(integer(3), x),
            mul(integer(-3), x)
        ),
        "0"
    );

    // ========================================================
    // Testes de Comutatividade das Operacoes
    // ========================================================

    std::cout
        << "============================================\n"
        << "22. Comutatividade das Operacoes\n"
        << "============================================\n\n";

    expect_equal(
        add(x, y),
        add(y, x),
        "Comutatividade da adição: x + y = y + x"
    );

    expect_equal(
        mul(x, y),
        mul(y, x),
        "Comutatividade da multiplicação: x * y = y * x"
    );

    expect_equal(
        mul(x, add(y, integer(2))),
        mul(add(y, integer(2)), x),
        "Comutatividade: x*(y+2) = (y+2)*x"
    );

    // ========================================================
    // Testes de Associatividade das Operacoes
    // ========================================================

    std::cout
        << "============================================\n"
        << "23. Associatividade das Operacoes\n"
        << "============================================\n\n";

    expect_equal(
        add(add(x, y), integer(3)),
        add(x, add(y, integer(3))),
        "Associatividade: (x+y)+3 = x+(y+3)"
    );

    expect_equal(
        mul(mul(x, y), integer(3)),
        mul(x, mul(y, integer(3))),
        "Associatividade: (x*y)*3 = x*(y*3)"
    );

    expect_equal(
        mul(x, add(y, integer(2))),
        mul(add(y, integer(2)), x),
        "Associatividade: x*(y+2) = (y+2)*x"
    );

    // ============================================
    // 24. Testes de Potências
    // ============================================

    std::cout << "\n============================================\n";
    std::cout << "24. Testes de Potências\n";
    std::cout << "============================================\n";

    // Teste 1: x^0 = 1
    expect_equal(
        pow(x, integer(0)),
        integer(1),
        "Potência: x^0 = 1"
    );

    // Teste 2: x^1 = x
    expect_equal(
        pow(x, integer(1)),
        x,
        "Potência: x^1 = x"
    );

    // Teste 3: (x^2)^3 = x^6
    auto x2 = pow(x, integer(2));
    auto x2_3 = pow(x2, integer(3));

    expect_equal(
        x2_3,
        pow(x, integer(6)),
        "Potência: (x^2)^3 = x^6"
    );

    // Teste 4: x*x^2 = x^3
    expect_equal(
        mul(x, x2),
        pow(x, integer(3)),
        "Potência: x*x^2 = x^3"
    );

// Teste 5: (x+y)^2 deve preservar a estrutura
auto xy2 = pow(add(x, y), integer(2));

std::cout << "[INFO] (x+y)^2 = " << str(xy2) << "\n";

    // ========================================================
    // Resumo
    // ========================================================

    std::cout
        << "\n============================================\n"
        << "                 SUMMARY\n"
        << "============================================\n\n";

    std::cout
        << "Testes executados : " << tests_run << "\n"
        << "Testes aprovados  : " << tests_passed << "\n"
        << "Testes falharam   : " << tests_failed << "\n\n";


    if (tests_failed == 0)
    {
        std::cout
            << "============================================\n"
            << "       TODOS OS TESTES PASSARAM!\n"
            << "============================================\n";

        return 0;
    }

    std::cout
        << "============================================\n"
        << "       ALGUNS TESTES FALHARAM!\n"
        << "============================================\n";

    return 1;
}