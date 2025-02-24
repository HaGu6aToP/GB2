#include <stdlib.h>
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
    if (fmpz_cmp(a, b) == 1){
        return a;
    } else {
        return b;
    }
}

void LCM(polynom res, fmpz** exp_f, fmpz** exp_g, slong nvars, fq_nmod_mpoly_ctx_t ctx){
    fq_nmod_mpoly_one(res, ctx);
    fmpz** new_exp = calloc(nvars, sizeof(fmpz*));
    for (int i = 0; i < nvars; i++){
        new_exp[i] = max(exp_g[i], exp_f[i]);
    }
    fq_nmod_mpoly_set_term_exp_fmpz(res, 0, new_exp, ctx);
    free(new_exp);
}

int S(polynom f, polynom g, fq_nmod_mpoly_ctx_t ctx){
    if (fq_nmod_mpoly_is_zero(f, ctx) && fq_nmod_mpoly_is_zero(g, ctx)){
        return 1;
    }
}

void buchberger_algoritm(basis new_basis, size_t new_basis_len, basis old_basis, size_t old_basis_len){

}

int main (void)
{
    // Инициализация простого поля
    fq_nmod_ctx_t field_ctx;
    fq_nmod_ctx_init_conway_ui(field_ctx, MEDIUM_PRIME, 1, "x");
    // fq_nmod_ctx_print(field_ctx); // Вывод информации о поле

    // Инициализация кольца многочленов над простым полем
    fq_nmod_mpoly_ctx_t poly_ring_ctx = {};
    // char** vars = calloc(3, sizeof(char[2]));
    // vars[0] = "x";
    // vars[1] = "y";
    // vars[2] = "z";

    const char* vars_array[] = {"x", "y", "z"};
    const char** vars = vars_array;

    slong nvars = 3; // Число переменных
    const ordering_t ordering = ORD_LEX; // Упорядочение
    fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, ordering, field_ctx);

    // Инициализация многочленов
    fq_nmod_mpoly_t p1, p2, p3, monom;
    fq_nmod_mpoly_init(p1, poly_ring_ctx);
    fq_nmod_mpoly_init(p2, poly_ring_ctx);
    fq_nmod_mpoly_init(p3, poly_ring_ctx);
    fq_nmod_mpoly_init(monom, poly_ring_ctx);

    // Задача 1
    fq_nmod_mpoly_set_str_pretty(p1, "x*y + x^2*z", vars, poly_ring_ctx);
    fq_nmod_mpoly_set_str_pretty(p2, "x*z + y*z^3", vars, poly_ring_ctx);
    fq_nmod_mpoly_set_str_pretty(p3, "y*z - y^2*z^3", vars, poly_ring_ctx);

    fq_nmod_mpoly_print_pretty(p1, vars, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(p2, vars, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(p3, vars, poly_ring_ctx);
    printf("\n");

    
    fmpz_t x, y, z, x2, y2, z2;
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(x2);
    fmpz_init(y2);
    fmpz_init(z2);

    fmpz *arr_exp[] = {x, y, z};
    fmpz** exp = arr_exp;

    fq_nmod_mpoly_get_term_exp_fmpz(exp, p3, 0, poly_ring_ctx);
    fmpz_print(exp[0]);
    fmpz_print(exp[1]);
    fmpz_print(exp[2]);

    fq_nmod_mpoly_one(p3, poly_ring_ctx);
    fq_nmod_mpoly_set_term_exp_fmpz(p3, 0, exp, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(p3, vars, poly_ring_ctx);
    printf("\n");

    fmpz *arr_exp_f[] = {x, y, z};
    fmpz *arr_exp_g[] = {x2, y2, z2};
    fmpz** exp_f = arr_exp_f;
    fmpz** exp_g = arr_exp_g;

    fq_nmod_mpoly_set_str_pretty(p1, "x^3*y^2 - x^2*y^3 + x", vars, poly_ring_ctx);
    fq_nmod_mpoly_set_str_pretty(p2, "2*x^4*y + y^2", vars, poly_ring_ctx);
    fq_nmod_mpoly_get_term_exp_fmpz(exp_f, p1, 0, poly_ring_ctx);
    fq_nmod_mpoly_get_term_exp_fmpz(exp_g, p2, 0, poly_ring_ctx);
    LCM(monom, exp_f, exp_g, nvars, poly_ring_ctx);
    fq_nmod_mpoly_print_pretty(monom, vars, poly_ring_ctx);
    printf("\n");

    // // Освобождение выделенных ресурсов
    fq_nmod_mpoly_ctx_clear(poly_ring_ctx);
    fq_nmod_ctx_clear(field_ctx);

    printf("\n");
    return 0;
}

