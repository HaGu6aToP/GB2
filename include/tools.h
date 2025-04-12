#pragma once
#include "headers.h"
#include "types.h" 

#define BUFFER_SIZE 1024

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