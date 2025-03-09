#include "headers.h"



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
#define VARS_COUNT 3

typedef fq_nmod_mpoly_ctx_t PolynomRing;
typedef fq_nmod_mpoly_t Polynom;
typedef fq_nmod_ctx_t Field;
typedef fmpz_t BIGINT;

ulong max(ulong a, ulong b){
    if (a > b) return a;
    else return b;
}

// Наименьшее общее кратное двух мономов
void LCM(Polynom monom, ulong nvars, const Polynom p1, const Polynom p2, const PolynomRing ctx){
    ulong exp_p1[nvars];
    ulong exp_p2[nvars];
    ulong exp_monom[nvars];

    fq_nmod_mpoly_get_term_exp_ui(exp_p1, p1, 0, ctx);
    fq_nmod_mpoly_get_term_exp_ui(exp_p2, p2, 0, ctx);

    for (int i = 0; i < nvars; i++)
        exp_monom[i] = max(exp_p1[i], exp_p2[i]);

    fq_nmod_mpoly_set_term_exp_ui(monom, 0, exp_monom, ctx);
}

// S полином
void S(Polynom S, ulong nvars, const Polynom p1, const Polynom p2, const PolynomRing ctx){
    fq_nmod_mpoly_t lcm;
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_one(lcm, ctx);
    LCM(lcm, nvars, p1, p2, ctx);

    fq_nmod_mpoly_t leading_monom_p1, leading_monom_p2, A;
    fq_nmod_mpoly_init(leading_monom_p1, ctx);
    fq_nmod_mpoly_init(leading_monom_p2, ctx);
    fq_nmod_mpoly_init(A, ctx);

    fq_nmod_mpoly_get_term(leading_monom_p1, p1, 0, ctx);
    fq_nmod_mpoly_get_term(leading_monom_p2, p2, 0, ctx);

    fq_nmod_mpoly_div(A, lcm, leading_monom_p1, ctx); // LCM(p1, p2) / LT(f)
    fq_nmod_mpoly_mul(leading_monom_p1, A, p1, ctx); // LCM(p1, p2) / LT(f) * p1

    fq_nmod_mpoly_div(A, lcm, leading_monom_p2, ctx); //LCM(p1, p2) / LT(f)
    fq_nmod_mpoly_mul(leading_monom_p2, A, p2, ctx); //LCM(p1, p2) / LT(f) * p2

    fq_nmod_mpoly_sub(S, leading_monom_p1, leading_monom_p2, ctx); // S

    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(leading_monom_p1, ctx);
    fq_nmod_mpoly_clear(leading_monom_p2, ctx);
    fq_nmod_mpoly_clear(A, ctx);
}

void log_S(Polynom S, ulong nvars, const Polynom p1, const Polynom p2, const PolynomRing ctx){
    fq_nmod_mpoly_t lcm;
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_one(lcm, ctx);
    LCM(lcm, nvars, p1, p2, ctx);

    printf("--------------------------\n");
    printf("LCM(f, g): ");
    fq_nmod_mpoly_print_pretty(lcm, NULL, ctx);
    printf("\n");
    
    fq_nmod_mpoly_t leading_monom_p1, leading_monom_p2, A;
    fq_nmod_mpoly_init(leading_monom_p1, ctx);
    fq_nmod_mpoly_init(leading_monom_p2, ctx);
    fq_nmod_mpoly_init(A, ctx);

    fq_nmod_mpoly_get_term(leading_monom_p1, p1, 0, ctx);
    fq_nmod_mpoly_get_term(leading_monom_p2, p2, 0, ctx);

    printf("LT(f): ");
    fq_nmod_mpoly_print_pretty(leading_monom_p1, NULL, ctx);
    printf("\n");

    printf("LT(g): ");
    fq_nmod_mpoly_print_pretty(leading_monom_p2, NULL, ctx);
    printf("\n");

    fq_nmod_mpoly_div(A, lcm, leading_monom_p1, ctx); // LCM(p1, p2) / LT(f)

    printf("LCM(p1, p2) / LT(f): ");
    fq_nmod_mpoly_print_pretty(A, NULL, ctx);
    printf("\n");

    fq_nmod_mpoly_mul(leading_monom_p1, A, p1, ctx); // LCM(p1, p2) / LT(f) * p1

    fq_nmod_mpoly_div(A, lcm, leading_monom_p2, ctx); //LCM(p1, p2) / LT(f)

    printf("LCM(p1, p2) / LT(g): ");
    fq_nmod_mpoly_print_pretty(A, NULL, ctx);
    printf("\n");

    fq_nmod_mpoly_mul(leading_monom_p2, A, p2, ctx); //LCM(p1, p2) / LT(f) * p2

    fq_nmod_mpoly_sub(S, leading_monom_p1, leading_monom_p2, ctx); // S

    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(leading_monom_p1, ctx);
    fq_nmod_mpoly_clear(leading_monom_p2, ctx);
    fq_nmod_mpoly_clear(A, ctx);
    printf("--------------------------\n");
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
    slong nvars = VARS_COUNT; // Число переменных
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
    const char* str_p2 = "3*x^4*y + y^2";//"x*z + y*z^3";
    const char* str_p3 = "x^3*y^2 - x^2*y^3 + x"; //"y*z - y^2*z^3";

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
    
    printf("\n");
    fq_nmod_mpoly_t s;
    fq_nmod_mpoly_init(s, poly_ring_ctx);
    S(s, VARS_COUNT, p3, p2, poly_ring_ctx);
    fq_nmod_mpoly_print_pretty(s, variables, poly_ring_ctx);
    printf("\n");

    printf("-----------GArray test------------------\n");
    GArray* basis = g_array_new(FALSE, FALSE, sizeof(fq_nmod_mpoly_struct));
    g_array_append_val(basis, p1);
    g_array_append_val(basis, p2);
    g_array_append_val(basis, p3);
    
    for (fq_nmod_mpoly_struct* i = (fq_nmod_mpoly_struct*)basis->data; i < (fq_nmod_mpoly_struct*)basis->data + basis->len; i++){
        fq_nmod_mpoly_print_pretty(i, variables, poly_ring_ctx);
        printf("\n");
    }
        
    printf("-----------------------------------------\n");

    g_array_remove_index(basis, 1);


    for (fq_nmod_mpoly_struct* i = (fq_nmod_mpoly_struct*)basis->data; i < (fq_nmod_mpoly_struct*)basis->data + basis->len; i++){
        fq_nmod_mpoly_print_pretty(i, variables, poly_ring_ctx);
        printf("\n");
    }

    fq_nmod_mpoly_struct* polynom = g_array_index(basis, fq_nmod_mpoly_t, 1);
    fq_nmod_mpoly_zero(polynom, poly_ring_ctx);

    printf("-----------------------------------------\n");
    for (fq_nmod_mpoly_struct* i = (fq_nmod_mpoly_struct*)basis->data; i < (fq_nmod_mpoly_struct*)basis->data + basis->len; i++){
        fq_nmod_mpoly_print_pretty(i, variables, poly_ring_ctx);
        printf("\n");
    }

    printf("---------cast to fq_nmd_mpoly_t----------n");
    fq_nmod_mpoly_t* G = (fq_nmod_mpoly_t*)basis->data;

    for (int i = 0; i < basis->len; i++){
        fq_nmod_mpoly_print_pretty(G[i], variables, poly_ring_ctx);
        printf("\n");
    }

    g_array_free(basis, TRUE);
    printf("-----------------------------------------\n");

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

