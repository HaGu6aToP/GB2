#include "headers.h"

#define SMALL_PRIME 7

ulong VARS_COUNT = 3;



ulong max(ulong a, ulong b){
    if (a > b) return a;
    else return b;
}

gint cmpPair(gconstpointer a, gconstpointer b){
    Pair *A = (Pair*)a, *B = (Pair*)b;
    // printf("\n---------------cmpPair---------------\n");
    // printf("A=(%ld, %ld), B=(%ld, %ld)\n", A->a, A->b, B->a, B->b);
    // printf("---------------------------------------\n");
    
    if (A->first < B->first) return -1;
    if (A->first > B->first) return 1;

    if (A->second < B->second) return -1;
    if (A->second > B->second) return 1;

    return 0;
}

// Наименьшее общее кратное ведущих мономов многочленов
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

    printf("Res: ");
    fq_nmod_mpoly_print_pretty(S, NULL, ctx);

    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(leading_monom_p1, ctx);
    fq_nmod_mpoly_clear(leading_monom_p2, ctx);
    fq_nmod_mpoly_clear(A, ctx);
    printf("\n--------------------------\n");
}

void log_B(GArray* B){
    printf("|B|=%d, B=\n", (int)B->len);
    for (Pair* i = (Pair*)B->data; i < (Pair*)B->data + B->len; i++){
        printf("(%ld, %ld)", i->first, i->second);
    }
    printf("\n\n");
}

void log_G(GArray* G, PolynomRing ctx){
    printf("\nG=\n");
    for (Polynom* p = (Polynom*)G->data; p < (Polynom*)G->data + G->len; p++){
        fq_nmod_mpoly_print_pretty(*p, NULL, ctx);
        printf("\n");
    }
    printf("\n");
}

// 1 - если критерий выполняется, 0 - не выполняется
int crit(GArray* G, GArray* B, ulong i, ulong j, const Polynom f, const Polynom g, const Polynom lcm, const PolynomRing ctx){
    int k;
    int res = 0;
    Pair pair_one, pair_two;
    fq_nmod_mpoly_t lt_h, Q, R;
    Polynom h;

    pair_one.first = i;
    pair_one.first = j;
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(Q, ctx);
    fq_nmod_mpoly_init(R, ctx);

    for(k = 0; k < G->len; k++){
        if (k == i || k == j) continue;
        h = g_array_index(G, Polynom, k);
        
        pair_one.second = k;
        pair_two.second = k;

        LT(lt_h, h, ctx);
        fq_nmod_mpoly_divrem(Q, R, lcm, lt_h, ctx);

        if ((g_array_binary_search(B, &pair_one, cmpPair, NULL) == FALSE) && (g_array_binary_search(B, &pair_two, cmpPair, NULL) == FALSE) && (fq_nmod_mpoly_is_zero(R, ctx))){
            res = 1;
            break;
        }
    }
    fq_nmod_mpoly_clear(lt_h, ctx);
    fq_nmod_mpoly_clear(Q, ctx);
    fq_nmod_mpoly_clear(R, ctx);
    return res;
}

Buchberger_result log_buchberger(const Basis basis, ulong t, const PolynomRing ctx){
    printf("\n------------------Buchberger-start-----------------\n");
    printf("<Инициализация>\n");
    
    ulong i, j, k, l;
    int flag1, flag2;
    flint_rand_t rand;
    Polynom f, g;
    fq_nmod_mpoly_t lcm, lt_f, lt_g, mul, S_polynom, S_mod_G;
    Pair pair;
    //----------------------------------------------------------------
    
    GArray* B = g_array_new(FALSE, FALSE, sizeof(Pair));
    GArray* G = g_array_new(FALSE, FALSE, sizeof(Polynom));
    flint_randinit(rand);
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_init(lt_f, ctx);
    fq_nmod_mpoly_init(lt_g, ctx);
    fq_nmod_mpoly_init(mul, ctx);
    fq_nmod_mpoly_init(S_polynom, ctx);
    fq_nmod_mpoly_init(S_mod_G, ctx);
    //----------------------------------------------------------------

    
    for(i = 0; i < t-1; i++){
        for(j = i + 1; j < t; j++){
            Pair pair = {i, j};
            g_array_append_val(B, pair);
        }
    }


    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(G, f);
    }

    log_B(B);
    
    log_G(G, ctx);


    printf("\n--------------------------Enter to Loop--------------------------\n");
    while(B->len != 0){
        k = n_randint(rand, B->len);
        pair = g_array_index(B, Pair, k);
        i = pair.first;
        j = pair.second;
        f = g_array_index(G, Polynom, i);
        g = g_array_index(G, Polynom, j);

        printf("Selected k: %ld\nCorresponding pair (i, j): {%ld, %ld}\n", k, i, j);

        printf("f: ");
        fq_nmod_mpoly_print_pretty(f, NULL, ctx);
        printf("\ng: ");
        fq_nmod_mpoly_print_pretty(g, NULL, ctx);
        printf("\n");

        LCM(lcm, f, g, ctx);
        printf("LCM(LT(f), LT(g)) = ");
        fq_nmod_mpoly_print_pretty(lcm, NULL, ctx);
        printf("\n");

        LT(lt_f, f, ctx);
        LT(lt_g, g, ctx);

        printf("LT(f) = ");
        fq_nmod_mpoly_print_pretty(lt_f, NULL, ctx);
        printf("\nLT(g) = ");
        fq_nmod_mpoly_print_pretty(lt_g, NULL, ctx);
        printf("\n");

        fq_nmod_mpoly_mul(mul, lt_f, lt_g, ctx);
        printf("LT(f)LT(g) = ");
        fq_nmod_mpoly_print_pretty(mul, NULL, ctx);
        printf("\n");

        flag1 = fq_nmod_mpoly_equal(lcm, mul, ctx);
        printf("crit call: (%ld, %ld)", i, j);
        flag2 = crit(G, B, i, j, f, g, lcm, ctx);

        printf("\nflag1=%d, flag2=%d\n", flag1, flag2);

        if ((flag1 == 0) && (flag2 == 0)){
            printf("Condition is met");
            printf("\n");

            S(S_polynom, f, g, ctx);
            printf("S(f, g) = ");
            fq_nmod_mpoly_print_pretty(S_polynom, NULL, ctx);
            printf("\n");

            Basis Q = init_empty_basis(G->len, ctx);
            Basis G_basis = from_garray(G);

            print_basis(G_basis, t, NULL, ctx);

            fq_nmod_mpoly_divrem_ideal(Q, S_mod_G, S_polynom, G_basis, t, ctx);
            // fq_nmod_mpoly_divrem_ideal((Polynom*)(Q->data), S_mod_G, S_polynom, (Polynom*)(G->data), G->len, ctx);
            printf("S(f, g) mod G = ");
            fq_nmod_mpoly_print_pretty(S_mod_G, NULL, ctx);
            printf("\n");
            
            free_basis(Q, t, ctx);
            flint_free(G_basis);

            if (fq_nmod_mpoly_is_zero(S_mod_G, ctx) == 0){
                printf("Adding new polynom to basis\n\n");
                Polynom new_polynom = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
                fq_nmod_mpoly_init(new_polynom, ctx);
                fq_nmod_mpoly_set(new_polynom, S_mod_G, ctx);
                // fq_nmod_mpoly_print_pretty(new_polynom, NULL, ctx);
                // printf("\n<--->\n");
                g_array_append_val(G, new_polynom);

                for(l = 0; l < t; l++){
                    Pair pair = {l, t};
                    g_array_append_val(B, pair);
                }

                t++;

                
            }
        }

        printf("Remove pair (%ld, %ld)\n\n", i, j);
        g_array_remove_index(B, k);
        log_B(B);
        printf("%d", G->len);
        log_G(G, ctx);
        // break;
    }


    //----------------------Освобождение ресурсов------------------------
    g_array_free(B, TRUE);
    flint_randclear(rand);
    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(lt_f, ctx);
    fq_nmod_mpoly_clear(lt_g, ctx);
    fq_nmod_mpoly_clear(mul, ctx);
    fq_nmod_mpoly_clear(S_polynom, ctx);
    fq_nmod_mpoly_clear(S_mod_G, ctx);
    printf("\n------------------Buchberger-end-----------------\n");
    
    Basis GBasis = from_garray(G);
    Buchberger_result res = {GBasis, G->len};
    g_array_free(G, TRUE);
    
    return res;
}

void reduce_groebner_basis(Basis basis, ulong len, PolynomRing ctx){

}

int is_groebner_basis(Basis basis, ulong len, PolynomRing ctx){
    ulong i, j;
    Polynom f, g;
    fq_nmod_mpoly_t S_polynom, S_mod_G;
    fq_nmod_mpoly_init(S_polynom, ctx);
    fq_nmod_mpoly_init(S_mod_G, ctx);
    int flag = 1;
    Basis Q = init_empty_basis(len, ctx);

    for(i = 0; i < len-1; i++){
        for(j = i; j < len; j++){
            f = basis[i];
            g = basis[j];
            log_S(S_polynom, f, g, ctx);
            fq_nmod_mpoly_divrem_ideal(Q, S_mod_G, S_polynom, basis, len, ctx);
            if (fq_nmod_mpoly_is_zero(S_mod_G, ctx) == 0){
                
                printf("\n\ni----------------------------is_groebner_basis----------------------------\n");
                printf("f: ");
                fq_nmod_mpoly_print_pretty(f, NULL, ctx);
                printf("\ng: ");
                fq_nmod_mpoly_print_pretty(g, NULL, ctx);
                printf("\nS(f, g): ");
                fq_nmod_mpoly_print_pretty(S_polynom, NULL, ctx);
                printf("\nS mod G: ");
                fq_nmod_mpoly_print_pretty(S_mod_G, NULL, ctx);
                printf("\n\n");

                free_basis(Q, len, ctx);
                fq_nmod_mpoly_clear(S_polynom, ctx);
                fq_nmod_mpoly_clear(S_mod_G, ctx);
                flag = 0;
                return flag;
            }
        }
    }

    fq_nmod_mpoly_clear(S_polynom, ctx);
    fq_nmod_mpoly_clear(S_mod_G, ctx);
    free_basis(Q, len, ctx);

    return flag;
}

Basis init_basis(ulong npolynoms, const char** strs, const char** vars, const PolynomRing ctx){
    Basis basis = flint_calloc(npolynoms, sizeof(Polynom));
    for (int i = 0; i < npolynoms; i++){
        basis[i] = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(basis[i], ctx);
        fq_nmod_mpoly_set_str_pretty(basis[i], strs[i], vars, ctx);
    }
    return basis;
}

Basis from_garray(GArray* g){
    Basis basis = flint_calloc(g->len, sizeof(Polynom));
    Polynom* p = (Polynom*)g->data;
    for (int i = 0; i < g->len; i++){
        basis[i] = *p;
        p += 1;
    }
    return basis;
}

Basis init_empty_basis(ulong npolynoms, const PolynomRing ctx){
    Basis basis = flint_calloc(npolynoms, sizeof(Polynom));
    for(int i = 0; i < npolynoms; i++){
        basis[i] = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(basis[i], ctx);
    }
    return basis;
}

void free_basis(Basis basis, ulong npolynoms, const PolynomRing ctx){
    for(int i = 0; i < npolynoms; i++){
        fq_nmod_mpoly_clear(basis[i], ctx);
        flint_free(basis[i]);
    }
    flint_free(basis);
}

void print_basis(const Basis basis, ulong npolnoms, const char** vars, const PolynomRing ctx){
    for(int i = 0; i < npolnoms; i++){
        fq_nmod_mpoly_print_pretty(basis[i], vars, ctx);
        printf("\n");
    }
    printf("\n");
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
    ulong npolynoms  = 3; // Число полиномов в системе
    ulong nvars = VARS_COUNT; // Число переменных
    const ordering_t order = ORD_DEGREVLEX; // Упорядочение
    fq_nmod_mpoly_ctx_t poly_ring_ctx = {};
    fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, order, field_ctx);

    const char* arr_variables[3] = {"x", "y", "z"};
    const char** variables = arr_variables;

    // char* str_p1 = "x*y + x^2*z";//"x^2 + y + z - 1";
    // char* str_p2 = "x*z + y*z^3";//"x + y^2 + z - 1";// "3*x^4*y + y^2";
    // char* str_p3 = "y*z - y^2*z^3";//"x + y + z^2 - 1";// "x^3*y^2 - x^2*y^3 + x"; 
    // const char** strs = flint_calloc(3, sizeof(char*));

    // strs[0] = str_p1;
    // strs[1] = str_p2;
    // strs[2] = str_p3;

    npolynoms = 2;
    char* str_p1 = "x^3 - 2*x*y";
    char* str_p2 ="x^2*y - 2*y^2 + x";
    const char** strs = flint_calloc(2, sizeof(char*));
    strs[0] = str_p1;
    strs[1] = str_p2;

    fq_nmod_mpoly_struct** basis = init_basis(npolynoms, strs, variables, poly_ring_ctx);
    print_basis(basis, npolynoms, variables, poly_ring_ctx);

    Buchberger_result GBasis = log_buchberger(basis, npolynoms, poly_ring_ctx);

    print_basis(GBasis.basis, GBasis.len, variables, poly_ring_ctx);

    int check = is_groebner_basis(GBasis.basis, GBasis.len, poly_ring_ctx);
    if (check == 1) printf("This is Groebner basis");
    else printf("This is not Groebner basis :c");

    free_basis(basis, npolynoms, poly_ring_ctx);
    free_basis(GBasis.basis, GBasis.len, poly_ring_ctx);

    // str_p1 = "x^3*y^2 - x^2*y^3 + x";
    // str_p2 = "3*x^4*y + y^2";
    // const char** test_strs = flint_calloc(2, sizeof(char*));
    // test_strs[0] = str_p1;
    // test_strs[1] = str_p2;
    // Basis test_basis = init_basis(2, test_strs, variables, poly_ring_ctx);
    // print_basis(test_basis, 2, variables, poly_ring_ctx);

    // fq_nmod_mpoly_t s;
    // fq_nmod_mpoly_init(s, poly_ring_ctx);
    // log_S(s, test_basis[0], test_basis[1], poly_ring_ctx);
    // fq_nmod_mpoly_print_pretty(s, variables, poly_ring_ctx);
    // printf("\n");

    // free_basis(test_basis, 2, poly_ring_ctx);

    
    
}