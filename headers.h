#include <stdlib.h>
#include <gmp.h>
#include "flint/flint.h"
#include "flint/fq_nmod.h" // Простое поле
#include "flint/fq_nmod_mpoly.h" // Кольцо многочленов нескольких переменных над конечным полем
#include "flint/fmpz.h" //
#include "flint/ulong_extras.h" // Для randint

// Динамический массив
#include "glib-2.0/glib.h"
#include "glib-2.0/glib/garray.h"

#define LT(res, f, ctx) fq_nmod_mpoly_get_term(res, f, 0, ctx)

typedef fq_nmod_mpoly_struct* Polynom;
typedef fq_nmod_mpoly_struct** Basis;
typedef fq_nmod_mpoly_ctx_struct* PolynomRing;

struct Pair{
    ulong first;
    ulong second;
};

struct Buchberger_result{
    Basis basis;
    ulong len;
};

typedef struct Pair Pair;
typedef struct Buchberger_result Buchberger_result;



Basis init_basis(ulong npolynoms, const char** strs, const char** vars, const PolynomRing ctx);
Basis init_empty_basis(ulong npolynoms, const PolynomRing ctx);
void free_basis(Basis basis, ulong npolynoms, const PolynomRing ctx);
void print_basis(const Basis basis, ulong npolnoms, const char** vars, const PolynomRing ctx);
Basis from_garray(GArray* g);
