#include "headers.h"


#define VARS_COUNT 3

#define SMALL_PRIME 7
#define MEDIUM_PRIME 43049
#define BIG_PRIME 56583254248148045179

typedef fq_nmod_mpoly_ctx_t PolynomRing;
typedef fq_nmod_mpoly_t Polynom;
typedef fq_nmod_ctx_t Field;

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
void LCM(Polynom monom, const Polynom p1, const Polynom p2, const PolynomRing ctx){
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
void S(Polynom S, const Polynom p1, const Polynom p2, const PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    fq_nmod_mpoly_t lcm;
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_one(lcm, ctx);
    LCM(lcm, p1, p2, ctx);

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

void log_S(Polynom S, const Polynom p1, const Polynom p2, const PolynomRing ctx){
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
        fq_nmod_mpoly_print_pretty(*p, NULL, ctx);
        printf("\n");
    }
}

int crit(GArray* G, GArray* B, ulong i, ulong j, Polynom f, Polynom g, PolynomRing ctx){
    ulong k;
    Polynom h;
    fq_nmod_mpoly_init(h, ctx);

    for (k = 0; k < B->len; k++){
        if (k == i || k == j)
            continue;
        fq_nmod_mpoly_set(h, g_array_index(G, Polynom, k), ctx);

    }
}

Polynom* buchberger(Polynom* basis, ulong t, PolynomRing ctx){
    ulong i, j, k, nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    GArray* B = g_array_new(FALSE, FALSE, sizeof(Pair));
    GArray* G = g_array_new(FALSE, FALSE, sizeof(Polynom));
    Pair curr_pair;
    Polynom f, g, mul, lcm, lt_f, lt_g;
    flint_rand_t random;

    _flint_rand_init_gmp(random);
    fq_nmod_mpoly_init(f, ctx);
    fq_nmod_mpoly_init(g, ctx);
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_init(mul, ctx);
    fq_nmod_mpoly_init(lt_f, ctx);
    fq_nmod_mpoly_init(lt_g, ctx);


    printf("\n-------------------------Buchberger-------------------------\n");

    for(i = 0; i < t-1; i++)
        for (j = i+1; j < t; j++){
            Pair pair = {i, j};
            g_array_append_val(B, pair);
        }

    log_B(B);
    Pair ss = {0, 2};
    int index = -1;
    int bool = g_array_binary_search(B, &ss, cmpPair, &index);
    printf("finded (%ld, %ld) index (%d) - %d\n", ss.a, ss.b, bool, index);
    printf("(%ld, %ld)\n", g_array_index(B, Pair, index).a, g_array_index(B, Pair, index).b);
    log_B(B);

    for (i = 0; i < nvars; i++){
        g_array_append_val(G, basis[i]);
    }

    log_G(G, ctx);

    printf("\n------------------------------------------------------------\n");

    while (B->len != 0){
        k = n_randint(random, B->len);
        curr_pair = g_array_index(B, Pair, k);
        printf("curr_pair = (%ld, %ld)\n", curr_pair.a, curr_pair.b);

        // *f = *g_array_index(G, Polynom, curr_pair.a);
        // *g = *g_array_index(G, Polynom, curr_pair.b);
        fq_nmod_mpoly_set(f, g_array_index(G, Polynom, curr_pair.a), ctx);
        fq_nmod_mpoly_set(g, g_array_index(G, Polynom, curr_pair.b), ctx);

        printf("f=");
        fq_nmod_mpoly_print_pretty(f, NULL, ctx);
        printf("\n");
        printf("g=");
        fq_nmod_mpoly_print_pretty(g, NULL, ctx);
        printf("\n");

        LCM(lcm, f, g, ctx);
        LT(lt_f, f, ctx);
        LT(lt_g, g, ctx);
        fq_nmod_mpoly_mul(mul, lt_f, lt_g, ctx);
        if (!fq_nmod_mpoly_equal(lcm, mul, ctx))
            1;

        break;
    }

    printf("\n------------------------------------------------------------\n");
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

    fq_nmod_mpoly_t* basis = flint_calloc(3, sizeof(fq_nmod_mpoly_t));
    fq_nmod_mpoly_init(basis[0], poly_ring_ctx);
    fq_nmod_mpoly_init(basis[1], poly_ring_ctx);
    fq_nmod_mpoly_init(basis[2], poly_ring_ctx);

    const char* arr_variables[3] = {"x", "y", "z"};
    const char** variables = arr_variables;

    const char* str_p1 = "x*y + x^2*z";
    const char* str_p2 = "3*x^4*y + y^2";//"x*z + y*z^3";
    const char* str_p3 = "x^3*y^2 - x^2*y^3 + x"; //"y*z - y^2*z^3";

    fq_nmod_mpoly_set_str_pretty(basis[0], str_p1, variables, poly_ring_ctx);
    fq_nmod_mpoly_set_str_pretty(basis[1], str_p2, variables, poly_ring_ctx); 
    fq_nmod_mpoly_set_str_pretty(basis[2], str_p3, variables, poly_ring_ctx);

    printf("Система:\n");
    fq_nmod_mpoly_print_pretty(basis[0], variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(basis[1], variables, poly_ring_ctx);
    printf("\n");
    fq_nmod_mpoly_print_pretty(basis[2], variables, poly_ring_ctx);
    printf("\n");

    buchberger(basis, 3, poly_ring_ctx);


}