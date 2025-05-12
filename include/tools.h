#pragma once
#include "headers.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "types.h" 

#define BUFFER_SIZE 1024

void free_variables(const char** variables, ulong nvars);
void get_variables(const char** variables, ulong nvars, const char* str);
void read_polinomials(Basis basis, ulong npoli,  const char** variables, PolynomRing ctx, FILE* file);
ulong max(ulong a, ulong b);
ulong max_poly_in_lst(const GArray* g, PolynomRing ctx);
ulong sum(ulong* arr, ulong len);
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
void poly_quick_sort(GArray* g, int first, int last, int rev, const PolynomRing ctx);
slong poly_binary_search(const GArray* g, const Polynom p, const PolynomRing ctx);
int is_poly_in_lst(const GArray* g, const Polynom p, const PolynomRing ctx);
void monom_lst_from_poly_lst(GArray* res, const GArray* g, const PolynomRing ctx);
void free_poly_lst(GArray* g, PolynomRing ctx);
void head_monom_lst_from_poly_lst(GArray* res, const GArray* g, const PolynomRing ctx);
void print_poly_lst(const GArray* lst, const PolynomRing ctx);

void* __calloc_poly_lst();
void* __calloc_poly();


#ifdef __cplusplus
}
#endif
