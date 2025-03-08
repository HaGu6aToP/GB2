#include <stdlib.h>
#include <gmp.h>
#include "flint/flint.h"
#include "flint/fq_nmod.h" // Простое поле
#include "flint/fq_nmod_mpoly.h" // Кольцо многочленов нескольких переменных над конечным полем
#include "flint/fmpz.h" //


#define SMALL_PRIME 7
#define MEDIUM_PRIME 43049
#define BIG_PRIME 56583254248148045179

// int main(int argc, char *argv[])
// {
//     ca_ctx_t ctx;
//     ca_t sqrt5, phi, psi, t, u;
//     fmpz_t n;

//     if (argc < 2)
//     {
//         flint_printf("usage: build/examples/binet [-limit B] n\n");
//         return 1;
//     }

//     fmpz_init(n);
//     fmpz_set_str(n, argv[argc-1], 10);

//     TIMEIT_ONCE_START
//     ca_ctx_init(ctx);

//     if (argc == 4)
//         ctx->options[CA_OPT_PREC_LIMIT] = atol(argv[2]);

//     ca_init(sqrt5, ctx);
//     ca_init(phi, ctx);
//     ca_init(psi, ctx);
//     ca_init(t, ctx);
//     ca_init(u, ctx);

//     ca_sqrt_ui(sqrt5, 5, ctx);

//     ca_add_ui(phi, sqrt5, 1, ctx);
//     ca_div_ui(phi, phi, 2, ctx);

//     ca_ui_sub(psi, 1, phi, ctx);

//     ca_pow_fmpz(t, phi, n, ctx);
//     ca_pow_fmpz(u, psi, n, ctx);

//     ca_sub(t, t, u, ctx);
//     ca_div(t, t, sqrt5, ctx);

//     ca_print(t, ctx);
//     flint_printf("\n");

//     ca_clear(sqrt5, ctx);
//     ca_clear(phi, ctx);
//     ca_clear(psi, ctx);
//     ca_clear(t, ctx);
//     ca_clear(u, ctx);
//     ca_ctx_clear(ctx);

//     fmpz_clear(n);
//     flint_printf("\n");
//     TIMEIT_ONCE_STOP
//     SHOW_MEMORY_USAGE

//     flint_cleanup();
//     return 0;
// }

typedef fq_nmod_mpoly_ctx_t* basis;
typedef fq_nmod_mpoly_t polynom;

fmpz* max(fmpz_t a, fmpz_t b){

}

// наименьшее общее кратное двух мономов
void LCM(polynom res, fmpz** exp_f, fmpz** exp_g, slong nvars, fq_nmod_mpoly_ctx_t ctx){

}

int S(polynom f, polynom g, fq_nmod_mpoly_ctx_t ctx){

}

void buchberger_algoritm(basis new_basis, size_t new_basis_len, basis old_basis, size_t old_basis_len){

}

int main (void)
{
    fmpz_t p; // порядок поля
    fmpz_init(p);
    fmpz_set_ui(p, SMALL_PRIME);

    printf("Порядок поля: ");
    fmpz_print(p);
    printf("\n");

    // Инициализация простого поля
    fq_nmod_ctx_t field_ctx;
    fq_nmod_ctx_init_conway(field_ctx, p, 1, "x");

    // Инициализация кольца многочленов над простым полем
    slong nvars = 3; // Число переменных
    const ordering_t order = ORD_LEX; // Упорядочение
    fq_nmod_mpoly_ctx_t poly_ring_ctx = {};
    fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, order, field_ctx);

    fq_nmod_mpoly_t p1, p2, p3;
    fq_nmod_mpoly_init(p1, poly_ring_ctx);
    fq_nmod_mpoly_init(p2, poly_ring_ctx);
    fq_nmod_mpoly_init(p3, poly_ring_ctx);

    const char* arr_variables[3] = {"x", "y", "z"};
    const char** variables = arr_variables;

    const char* str_p1 = "x*y + x^2*z";
    const char* str_p2 = "x*z + y*z^3";
    const char* str_p3 = "y*z - y^2*z^3";

    fq_nmod_mpoly_set_str_pretty(p1, str_p1, variables, poly_ring_ctx);
    fq_nmod_mpoly_set_str_pretty(p2, str_p2, variables, poly_ring_ctx); 
    fq_nmod_mpoly_set_str_pretty(p3, str_p3, variables, poly_ring_ctx);

    printf("Система:\n");
    fq_nmod_mpoly_print_pretty(p1, variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(p2, variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(p3, variables, poly_ring_ctx);
    printf("\n");

    // Освобождение ресурсов
    fq_nmod_mpoly_clear(p1, poly_ring_ctx);
    fq_nmod_mpoly_clear(p2, poly_ring_ctx);
    fq_nmod_mpoly_clear(p3, poly_ring_ctx);
    fq_nmod_mpoly_ctx_clear(poly_ring_ctx);
    fq_nmod_ctx_clear(field_ctx);



    flint_cleanup();
    printf("\n");
    return 0;
}

