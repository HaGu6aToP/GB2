#pragma once

#include <stdlib.h>
#include <gmp.h>
#include "flint/flint.h"
#include "flint/fq_nmod.h" // Простое поле
#include "flint/fq_nmod_mpoly.h" // Кольцо многочленов нескольких переменных над конечным полем
#include "flint/fmpz.h" //
#include "flint/ulong_extras.h" // Для randint
#include "flint/thread_support.h"

// Динамический массив
#include "glib-2.0/glib.h"
#include "glib-2.0/glib/garray.h"

#define LT(res, f, ctx) fq_nmod_mpoly_get_term(res, f, 0, ctx)
#define BUFFER_SIZE 1024

typedef fq_nmod_mpoly_struct* Polynom;
typedef fq_nmod_mpoly_struct** Basis;
typedef fq_nmod_mpoly_ctx_struct* PolynomRing;

struct Pair{
    ulong first;
    ulong second;
};

struct SPair{
    Polynom poly;
    ulong first;
    ulong second;
};

struct Buchberger_result{
    Basis basis;
    ulong len;
};

struct thread_buchberger_data_t{
    ulong* finished_threads;
    Polynom S_poly;
    Basis F;
    ulong npoli;
    PolynomRing ctx;
    int completed;
};

typedef struct Pair Pair;
typedef struct SPair SPair;
typedef struct Buchberger_result Buchberger_result;
typedef struct thread_buchberger_data_t thread_buchberger_data_t;

extern ulong NO_OF_IRRED;


// basis_tools
Basis init_basis(ulong npoli, const char** strs, const char** vars, const PolynomRing ctx);
Basis init_empty_basis(ulong npoli, const PolynomRing ctx);
void free_basis(Basis basis, ulong npoli, const PolynomRing ctx);
void print_basis(const Basis basis, ulong npoli, const char** vars, const PolynomRing ctx);
Basis from_garray(GArray* g);

// tools
void free_variables(const char** variables, ulong nvars);
void get_variables(const char** variables, ulong nvars, const char* str);
void read_polinomials(Basis basis, ulong npoli,  const char** variables, PolynomRing ctx, FILE* file);
ulong max(ulong a, ulong b);
gint cmpPair(gconstpointer a, gconstpointer b);
void log_B(GArray* B);
void log_G(GArray* G, PolynomRing ctx);
int parseInt(char* chars);
int powInt(int x, int y);
void print_poly(const char* header, Polynom p, const char** vars, PolynomRing ctx);
void init_SPair(SPair* pspair, PolynomRing ctx);
void set_SPair(SPair* pspair, Polynom p, ulong first, ulong second, PolynomRing ctx);
void free_SPair(SPair* pspair, PolynomRing ctx);
void copy_SPair(SPair* pspair, const SPair* resourse);
void print_ulong_garray(GArray* g);
void quick_sort(ulong* arr, int left, int right);

// buchberger
void LCM(Polynom monom, const Polynom p1, const Polynom p2, const PolynomRing ctx);
void S(Polynom S, const Polynom p1, const Polynom p2, const PolynomRing ctx);
void log_S(Polynom S, const Polynom p1, const Polynom p2, const PolynomRing ctx);
int crit(GArray* G, GArray* B, ulong i, ulong j, const Polynom f, const Polynom g, const Polynom lcm, const PolynomRing ctx);
Buchberger_result log_buchberger(const Basis basis, ulong t, const PolynomRing ctx);
Buchberger_result buchberger(const Basis basis, ulong t, const PolynomRing ctx);
int is_groebner_basis(Basis basis, ulong len, PolynomRing ctx);
Buchberger_result buchberger_v2(const Basis basis, ulong t, const PolynomRing ctx);
Buchberger_result log_buchberger_v2(const Basis basis, ulong t, const PolynomRing ctx);
Buchberger_result buchberger_v2_1(const Basis basis, ulong t, const PolynomRing ctx);
int find_min(GArray* P, PolynomRing ctx);
SPair find_min_v1(GArray* F, GArray* P, PolynomRing ctx);
int log_find_min(GArray* P, PolynomRing ctx);

Buchberger_result threaded_buchberger(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx);
Buchberger_result log_threaded_buchberger(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx);
void my_exit(thread_buchberger_data_t* data);