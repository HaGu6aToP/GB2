#pragma once
#include "headers.h"
#ifdef __cplusplus
extern "C" {
#endif

#include "types.h"

Basis init_basis(ulong npoli, const char** strs, const char** vars, const PolynomRing ctx);
Basis init_empty_basis(ulong npoli, const PolynomRing ctx);
void free_basis(Basis basis, ulong npoli, const PolynomRing ctx);
void print_basis(const Basis basis, ulong npoli, const char** vars, const PolynomRing ctx);
Basis from_garray(GArray* g);

#ifdef __cplusplus
}
#endif