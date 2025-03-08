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

struct Pair{
    ulong a;
    ulong b;
};

typedef struct Pair Pair;