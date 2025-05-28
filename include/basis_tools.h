#pragma once
#include "headers.h"
#ifdef __cplusplus
extern "C" {
#endif

#include "types.h"

// Выделяет память для массива из npoli многочленов и инициализирует их из строки strs с переменными vars
Basis init_basis(ulong npoli, const char** strs, const char** vars, const PolynomRing ctx);

// Выделить память для массива из npolu многочленов и инициализирует их нулевыми многочленами
Basis init_empty_basis(ulong npoli, const PolynomRing ctx);

// Освобождение памяти выделенной для базиса
void free_basis(Basis basis, ulong npoli, const PolynomRing ctx);

// Вывести базис в консоль
void print_basis(const Basis basis, ulong npoli, const char** vars, const PolynomRing ctx);

Basis from_garray(GArray* g);

#ifdef __cplusplus
}
#endif
