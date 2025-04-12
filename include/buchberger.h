#pragma once

#include "headers.h"
#include "types.h"
#include "flint/fmpz.h" //
#include "flint/ulong_extras.h" // Для randint
#include "flint/thread_support.h"

#define LT(res, f, ctx) fq_nmod_mpoly_get_term(res, f, 0, ctx)


extern ulong NO_OF_IRRED;

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
Buchberger_result threaded_buchberger_v2(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx);
Buchberger_result log_threaded_buchberger(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx);
void my_exit(thread_buchberger_data_t* data);