#include "headers.h"

#define SMALL_PRIME 7

ulong VARS_COUNT = 3;

typedef fq_nmod_mpoly_struct Polynom;
typedef fq_nmod_mpoly_struct* Basis;
typedef fq_nmod_mpoly_ctx_t PolynomRing;

ulong max(ulong a, ulong b){
    if (a > b) return a;
    else return b;
}

gint cmpPair(gconstpointer a, gconstpointer b){
    Pair *A = (Pair*)a, *B = (Pair*)b;
    // printf("\n---------------cmpPair---------------\n");
    // printf("A=(%ld, %ld), B=(%ld, %ld)\n", A->a, A->b, B->a, B->b);
    // printf("---------------------------------------\n");
    
    if (A->a < B->a) return -1;
    if (A->a > B->a) return 1;

    if (A->b < B->b) return -1;
    if (A->b > B->b) return 1;

    return 0;
}

// Наименьшее общее кратное двух мономов
void LCM(Polynom* monom, const Polynom* p1, const Polynom* p2, const PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    ulong exp_p1[nvars];
    ulong exp_p2[nvars];
    ulong exp_monom[nvars];

    fq_nmod_mpoly_get_term_exp_ui(exp_p1, p1, 0, ctx);
    fq_nmod_mpoly_get_term_exp_ui(exp_p2, p2, 0, ctx);

    for (int i = 0; i < nvars; i++)
        exp_monom[i] = max(exp_p1[i], exp_p2[i]);

    fq_nmod_mpoly_one(monom, ctx);
    fq_nmod_mpoly_set_term_exp_ui(monom, 0, exp_monom, ctx);
}

// S полином
void S(Polynom* S, const Polynom* p1, const Polynom* p2, const PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    Polynom lcm;
    fq_nmod_mpoly_init(&lcm, ctx);
    fq_nmod_mpoly_one(&lcm, ctx);
    LCM(&lcm, p1, p2, ctx);

    Polynom leading_monom_p1, leading_monom_p2, A;
    fq_nmod_mpoly_init(&leading_monom_p1, ctx);
    fq_nmod_mpoly_init(&leading_monom_p2, ctx);
    fq_nmod_mpoly_init(&A, ctx);

    fq_nmod_mpoly_get_term(&leading_monom_p1, p1, 0, ctx);
    fq_nmod_mpoly_get_term(&leading_monom_p2, p2, 0, ctx);

    fq_nmod_mpoly_div(&A, &lcm, &leading_monom_p1, ctx); // LCM(p1, p2) / LT(f)
    fq_nmod_mpoly_mul(&leading_monom_p1, &A, p1, ctx); // LCM(p1, p2) / LT(f) * p1

    fq_nmod_mpoly_div(&A, &lcm, &leading_monom_p2, ctx); //LCM(p1, p2) / LT(f)
    fq_nmod_mpoly_mul(&leading_monom_p2, &A, p2, ctx); //LCM(p1, p2) / LT(f) * p2

    fq_nmod_mpoly_sub(S, &leading_monom_p1, &leading_monom_p2, ctx); // S

    fq_nmod_mpoly_clear(&lcm, ctx);
    fq_nmod_mpoly_clear(&leading_monom_p1, ctx);
    fq_nmod_mpoly_clear(&leading_monom_p2, ctx);
    fq_nmod_mpoly_clear(&A, ctx);
}

void log_S(Polynom* S, const Polynom* p1, const Polynom* p2, const PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);

    fq_nmod_mpoly_t lcm;
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_one(lcm, ctx);
    LCM(lcm, p1, p2, ctx);

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

void log_B(GArray* B){
    printf("%d\nB=\n", (int)B->len);
    for (Pair* i = (Pair*)B->data; i < (Pair*)B->data + B->len; i++){
        printf("(%ld, %ld)", i->a, i->b);
    }
    printf("\n\n");
}

void log_G(GArray* G, PolynomRing ctx){
    printf("\nG=\n");
    for (Polynom* p = (Polynom*)G->data; p < (Polynom*)G->data + G->len; p++){
        fq_nmod_mpoly_print_pretty(p, NULL, ctx);
        printf("\n");
    }
}

int crit(GArray* G, GArray* B, ulong i, ulong j, const Polynom* f, const Polynom* g, const PolynomRing ctx){
    ulong k;
    int res = 0;
    Polynom h, lcm, lt_h, div;
    Pair pair_i, pair_j;
    pair_i.a = i;
    pair_j.a = j;
    fq_nmod_mpoly_init(&h, ctx);
    fq_nmod_mpoly_init(&lcm, ctx);
    fq_nmod_mpoly_init(&lt_h, ctx);
    fq_nmod_mpoly_init(&div, ctx);

    for (k = 0; k < B->len; k++){
        if (k == i || k == j)
            continue;
        
        pair_i.b = k;
        pair_j.b = k;
        fq_nmod_mpoly_set(&h, &g_array_index(G, Polynom, k), ctx);
        LT(&lt_h, &h, ctx);
        LCM(&lcm, f, g, ctx);
        

        if (g_array_binary_search(B, &pair_i, cmpPair, NULL) == 0 && g_array_binary_search(B, &pair_j, cmpPair, NULL) == 0 && fq_nmod_mpoly_divides(&div, &lcm, &lt_h, ctx) == 1){
            res = 1;
            break;
        }
    }

    fq_nmod_mpoly_clear(&h, ctx);
    fq_nmod_mpoly_clear(&lcm, ctx);
    fq_nmod_mpoly_clear(&lt_h, ctx);
    fq_nmod_mpoly_clear(&div, ctx);
    // printf("res = %d\n", res);
    return res;
}

Basis log_buchberger(Basis basis, ulong t, PolynomRing ctx){
    printf("------------------Buchberger------------------");
    printf("<Инициализация>");
}

void main(void){
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

    fq_nmod_mpoly_struct* basis = flint_calloc(3, sizeof(fq_nmod_mpoly_struct));
    fq_nmod_mpoly_init(&basis[0], poly_ring_ctx);
    fq_nmod_mpoly_init(&basis[1], poly_ring_ctx);
    fq_nmod_mpoly_init(&basis[2], poly_ring_ctx);

    const char* arr_variables[3] = {"x", "y", "z"};
    const char** variables = arr_variables;

    const char* str_p1 = "x*y + x^2*z";
    const char* str_p2 = "x*z + y*z^3";// "3*x^4*y + y^2";
    const char* str_p3 = "y*z - y^2*z^3";// "x^3*y^2 - x^2*y^3 + x"; 

    fq_nmod_mpoly_set_str_pretty(&basis[0], str_p1, variables, poly_ring_ctx);
    fq_nmod_mpoly_set_str_pretty(&basis[1], str_p2, variables, poly_ring_ctx); 
    fq_nmod_mpoly_set_str_pretty(&basis[2], str_p3, variables, poly_ring_ctx);

    printf("Система:\n");
    fq_nmod_mpoly_print_pretty(&basis[0], variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(&basis[1], variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(&basis[2], variables, poly_ring_ctx);
    printf("\n");
    
    // log_buchberger(basis, 3, poly_ring_ctx);

    fq_nmod_mpoly_struct* Q = flint_calloc(3, sizeof(fq_nmod_mpoly_struct));
    fq_nmod_mpoly_t s, Sp;
    fq_nmod_mpoly_init(s, poly_ring_ctx);
    fq_nmod_mpoly_init(Sp, poly_ring_ctx);
    S(Sp, &basis[0], &basis[1], poly_ring_ctx);

    printf("Система:\n");
    fq_nmod_mpoly_print_pretty(&basis[0], variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(&basis[1], variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(&basis[2], variables, poly_ring_ctx);
    printf("\n");

    fq_nmod_mpoly_divrem_ideal(&Q, s, Sp, &basis, 3, poly_ring_ctx);

    fq_nmod_mpoly_print_pretty(s, variables, poly_ring_ctx);
    
    for (int i = 0; i < 3; i++)
        fq_nmod_mpoly_init(&Q[i], poly_ring_ctx);

    fq_nmod_mpoly_clear(&basis[0], poly_ring_ctx);
    fq_nmod_mpoly_clear(&basis[1], poly_ring_ctx);
    fq_nmod_mpoly_clear(&basis[2], poly_ring_ctx);

    for (int i = 0; i < 3; i++)
        fq_nmod_mpoly_clear(&Q[i], poly_ring_ctx);


}