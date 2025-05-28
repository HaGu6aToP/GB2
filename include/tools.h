#pragma once
#include "headers.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "types.h" 

#define BUFFER_SIZE 1024

// Освободить память выделенную для variables
void free_variables(const char** variables, ulong nvars);

// Получить переменные из строки. Например "x y z"
void get_variables(const char** variables, ulong nvars, const char* str);

// Прочитать систему многочленов из файла
void read_polinomials(Basis basis, ulong npoli,  const char** variables, PolynomRing ctx, FILE* file);

// Максимум из двух чисел
ulong max(ulong a, ulong b);

// Максимальный полином в массиве
ulong max_poly_in_lst(const GArray* g, PolynomRing ctx);

// Сумма элементов масисва
ulong sum(ulong* arr, ulong len);

// Отношение сравнение S-пар для сортировки
gint cmpPair(gconstpointer a, gconstpointer b);

void log_B(GArray* B);
void log_G(GArray* G, PolynomRing ctx);

// int из строки
int parseInt(char* chars);
int powInt(int x, int y);

// Вывести многочлен в консоль
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
