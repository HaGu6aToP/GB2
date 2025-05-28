#include "buchberger.h"
#include "tools.h"
#include "basis_tools.h"

#include <pthread.h>
#include <unistd.h>

ulong NO_OF_IRRED = 1;


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


int crit(GArray* G, GArray* B, ulong i, ulong j, const Polynom f, const Polynom g, const Polynom lcm, const PolynomRing ctx){
    ulong k;
    int flag1, flag2, flag3;
    int res = 0;
    Pair pair_one, pair_two;
    fq_nmod_mpoly_t lt_h, Q, R;
    Polynom h;

    // pair_one.first = i;
    // pair_one.first = j;
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(Q, ctx);
    fq_nmod_mpoly_init(R, ctx);

    // log_B(B);
    // log_G(G, ctx);

    for(k = 0; k < G->len; k++){
        if (k == i || k == j) continue;
        h = g_array_index(G, Polynom, k);

        // printf("\nSelected k: %ld", k);
        // printf("\nh=");
        // fq_nmod_mpoly_print_pretty(h, NULL, ctx);
        // printf("\n");
        
        
        if (i < k){
            pair_one.first = i;
            pair_one.second = k;
        } else {
            pair_one.first = k;
            pair_one.second = i;
        }

        if (j < k){
            pair_two.first = j;
            pair_two.second = k;
        } else {
            pair_two.first = k;
            pair_two.second = j;
        }

        LT(lt_h, h, ctx);
        fq_nmod_mpoly_divrem(Q, R, lcm, lt_h, ctx);

        g_array_sort(B, cmpPair);
        flag1 = g_array_binary_search(B, &pair_one, cmpPair, NULL);
        flag2 = g_array_binary_search(B, &pair_two, cmpPair, NULL);
        flag3 = fq_nmod_mpoly_is_zero(R, ctx);

        // printf("Найдена ли пара (%ld, %ld): %d\n", i, k, flag1);
        // printf("Найдена ли пара (%ld, %ld): %d\n", j, k, flag2);
        // printf("Остаток от деления: ");
        // fq_nmod_mpoly_print_pretty(R, NULL, ctx);
        // printf("\n");

        if ((flag1 == 0) && ( flag2 == 0) && (fq_nmod_mpoly_is_zero(R, ctx))){
            res = 1;
            break;
        }
    }
    fq_nmod_mpoly_clear(lt_h, ctx);
    fq_nmod_mpoly_clear(Q, ctx);
    fq_nmod_mpoly_clear(R, ctx);
    return res;
}

int log_crit(GArray* G, GArray* B, ulong i, ulong j, const Polynom f, const Polynom g, const Polynom lcm, const PolynomRing ctx){
    ulong k;
    int flag1, flag2, flag3;
    int res = 0;
    Pair pair_one, pair_two;
    fq_nmod_mpoly_t lt_h, Q, R;
    Polynom h;

    // pair_one.first = i;
    // pair_one.first = j;
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(Q, ctx);
    fq_nmod_mpoly_init(R, ctx);

    log_B(B);
    log_G(G, ctx);

    for(k = 0; k < G->len; k++){
        if (k == i || k == j) continue;
        h = g_array_index(G, Polynom, k);

        printf("\nSelected k: %ld", k);
        printf("\nh=");
        fq_nmod_mpoly_print_pretty(h, NULL, ctx);
        printf("\n");
        
        
        if (i < k){
            pair_one.first = i;
            pair_one.second = k;
        } else {
            pair_one.first = k;
            pair_one.second = i;
        }

        if (j < k){
            pair_two.first = j;
            pair_two.second = k;
        } else {
            pair_two.first = k;
            pair_two.second = j;
        }

        LT(lt_h, h, ctx);
        fq_nmod_mpoly_divrem(Q, R, lcm, lt_h, ctx);

        g_array_sort(B, cmpPair);
        flag1 = g_array_binary_search(B, &pair_one, cmpPair, NULL);
        flag2 = g_array_binary_search(B, &pair_two, cmpPair, NULL);
        flag3 = fq_nmod_mpoly_is_zero(R, ctx);

        printf("Найдена ли пара (%ld, %ld): %d\n", i, k, flag1);
        printf("Найдена ли пара (%ld, %ld): %d\n", j, k, flag2);
        printf("Остаток от деления: ");
        fq_nmod_mpoly_print_pretty(R, NULL, ctx);
        printf("\n");

        if ((flag1 == 0) && ( flag2 == 0) && (fq_nmod_mpoly_is_zero(R, ctx))){
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
        printf("\n------------crit call------------\n");
        log_G(G, ctx);
        log_B(B);
        printf("i=%ld, j=%ld\nf=", i, j);
        fq_nmod_mpoly_print_pretty(f, NULL, ctx);
        printf("\ng=");
        fq_nmod_mpoly_print_pretty(g, NULL, ctx);
        printf("\nLCM(f, g)=");
        fq_nmod_mpoly_print_pretty(lcm, NULL, ctx);
        printf("\n");
        flag2 = log_crit(G, B, i, j, f, g, lcm, ctx);

        printf("\n---------------------------------\n");

        printf("\nflag1=%d, flag2=%d\n", flag1, flag2);

        if ((flag1 == 0) && (flag2 == 0)){
            printf("Condition is met");
            printf("\nf=");
            fq_nmod_mpoly_print_pretty(f, NULL, ctx);
            printf("\ng=");
            fq_nmod_mpoly_print_pretty(g, NULL, ctx);
            printf("\nLCM(f, g)=");
            fq_nmod_mpoly_print_pretty(lcm, NULL, ctx);
            printf("\n");

            S(S_polynom, f, g, ctx);
            printf("S(f, g) = ");
            fq_nmod_mpoly_print_pretty(S_polynom, NULL, ctx);
            printf("\n");

            Basis Q = init_empty_basis(G->len, ctx);
            Basis G_basis = from_garray(G);

            // print_basis(G_basis, t, NULL, ctx);

            fq_nmod_mpoly_divrem_ideal(Q, S_mod_G, S_polynom, G_basis, G->len, ctx);
            // fq_nmod_mpoly_divrem_ideal((Polynom*)(Q->data), S_mod_G, S_polynom, (Polynom*)(G->data), G->len, ctx);
            printf("S(f, g) mod G = ");
            fq_nmod_mpoly_print_pretty(S_mod_G, NULL, ctx);
            printf("\n");
            
            free_basis(Q, G->len, ctx);
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
            S(S_polynom, f, g, ctx);
            fq_nmod_mpoly_divrem_ideal(Q, S_mod_G, S_polynom, basis, len, ctx);
            if (fq_nmod_mpoly_is_zero(S_mod_G, ctx) == 0){
                
                // printf("\n\ni----------------------------is_groebner_basis----------------------------\n");
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



Buchberger_result buchberger(const Basis basis, ulong t, const PolynomRing ctx){
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

    while(B->len != 0){
        // printf("B->len=%d\n", B->len);

        k = n_randint(rand, B->len);
        pair = g_array_index(B, Pair, k);
        i = pair.first;
        j = pair.second;
        f = g_array_index(G, Polynom, i);
        g = g_array_index(G, Polynom, j);

        LCM(lcm, f, g, ctx);
        LT(lt_f, f, ctx);
        LT(lt_g, g, ctx);

        fq_nmod_mpoly_mul(mul, lt_f, lt_g, ctx);

        flag1 = fq_nmod_mpoly_equal(lcm, mul, ctx);
        flag2 = crit(G, B, i, j, f, g, lcm, ctx);

        if ((flag1 == 0) && (flag2 == 0)){
            S(S_polynom, f, g, ctx);

            Basis Q = init_empty_basis(G->len, ctx);
            Basis G_basis = from_garray(G);

            // G_basis = (Basis)G->data;

            fq_nmod_mpoly_divrem_ideal(Q, S_mod_G, S_polynom, G_basis, G->len, ctx);

            
            free_basis(Q, G->len, ctx);
            flint_free(G_basis);

            if (fq_nmod_mpoly_is_zero(S_mod_G, ctx) == 0){
                Polynom new_polynom = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
                fq_nmod_mpoly_init(new_polynom, ctx);
                fq_nmod_mpoly_set(new_polynom, S_mod_G, ctx);
                g_array_append_val(G, new_polynom);

                for(l = 0; l < t; l++){
                    Pair pair = {l, t};
                    g_array_append_val(B, pair);
                }

                t++;

                
            }
        }

        g_array_remove_index(B, k);
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
    
    Basis GBasis = from_garray(G);
    Buchberger_result res = {GBasis, G->len};
    g_array_free(G, TRUE);
    
    return res;
}


void get_data_find_min_v1(Pair* ppair, GArray* F, int* i, int* j, Polynom f, Polynom g, Polynom S_poly, PolynomRing ctx){
    *i = ppair->first;
    *j = ppair->second;
    f = g_array_index(F, Polynom, *i);
    g = g_array_index(F, Polynom, *j);
    S(S_poly, f, g, ctx);
}

int log_find_min(GArray* P, PolynomRing ctx){
    fq_nmod_mpoly_t lt_min, lt_S;
    int res = 0;
    SPair* pspair;
    int i;
    printf("------------------find_min------------------\n");
//----------------------------------------------------
    fq_nmod_mpoly_init(lt_min, ctx);
    fq_nmod_mpoly_init(lt_S, ctx);
//----------------------------------------------------
    if (P->len > 0){
        pspair = (SPair*)P->data;
        LT(lt_min, pspair->poly, ctx);
        res = 0;
        pspair++;

        for(i = 1; i < P->len; i++){
            print_poly("poly:", pspair->poly, NULL, ctx);
            LT(lt_S, pspair->poly, ctx);
            printf("i=%d, %d\n", i, P->len);
            print_poly("curr min:", lt_min, NULL, ctx);
            print_poly("curr poly:", lt_S, NULL, ctx);
            printf("cmp=%d\n", fq_nmod_mpoly_cmp(lt_min, lt_S, ctx));
            if (fq_nmod_mpoly_cmp(lt_min, lt_S, ctx) == 1){
                fq_nmod_mpoly_set(lt_min, lt_S, ctx);
                res = i;
            }
            pspair++;
        }
    }
//----------------------------------------------------
    fq_nmod_mpoly_clear(lt_min, ctx);
    fq_nmod_mpoly_clear(lt_S, ctx);
    printf("--------------------------------------------\n");
//----------------------------------------------------
    return res;
}

int find_min(GArray* P, PolynomRing ctx){
    fq_nmod_mpoly_t lt_min, lt_S;
    int res = 0;
    SPair* pspair;
    int i;
//----------------------------------------------------
    fq_nmod_mpoly_init(lt_min, ctx);
    fq_nmod_mpoly_init(lt_S, ctx);
//----------------------------------------------------
    if (P->len > 0){
        pspair = (SPair*)P->data;
        LT(lt_min, pspair->poly, ctx);
        res = 0;
        pspair++;

        for(i = 1; i < P->len; i++){
            LT(lt_S, pspair->poly, ctx);
            if (fq_nmod_mpoly_cmp(lt_min, lt_S, ctx) == 1){
                fq_nmod_mpoly_set(lt_min, lt_S, ctx);
                res = i;
            }
            pspair++;
        }
    }
//----------------------------------------------------
    fq_nmod_mpoly_clear(lt_min, ctx);
    fq_nmod_mpoly_clear(lt_S, ctx);
//----------------------------------------------------
    return res;
}


SPair find_min_v1(GArray* F, GArray* P, PolynomRing ctx){
    fq_nmod_mpoly_t S_poly, lt_min, lt_S;
    SPair res;
    Polynom f, g;
    Pair* ppair;
    int i, j, k;
//----------------------------------------------------
    fq_nmod_mpoly_init(S_poly, ctx);
    fq_nmod_mpoly_init(lt_min, ctx);
    fq_nmod_mpoly_init(lt_S, ctx);
//----------------------------------------------------
    ppair = (Pair*)P->data;
    get_data_find_min_v1(ppair, F, &i, &j, f, g, S_poly, ctx);
    if (fq_nmod_mpoly_is_zero(S_poly, ctx) == 0) LT(lt_min, S_poly, ctx);

    // res.poly = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    init_SPair(&res, ctx);
    set_SPair(&res, S_poly, i, j, ctx);
    ppair++;
    
    for(k = 1; k < P->len; k++){
        // i = ppair->first;
        // j = ppair->second;
        // f = g_array_index(F, Polynom, i);
        // g = g_array_index(F, Polynom, j);
        // S(S_poly, f, g, ctx);
        get_data_find_min_v1(ppair, F, &i, &j, f, g, S_poly, ctx);


        if (fq_nmod_mpoly_is_zero(S_poly, ctx) == 0) LT(lt_S, S_poly, ctx);

        if (fq_nmod_mpoly_cmp(lt_min, lt_S, ctx) == 1){
            fq_nmod_mpoly_set(lt_min, lt_S, ctx);
            set_SPair(&res, S_poly, i, j, ctx);
        }

        ppair++;
    }
//----------------------------------------------------
    fq_nmod_mpoly_clear(S_poly, ctx);
    fq_nmod_mpoly_clear(lt_min, ctx);
    fq_nmod_mpoly_clear(lt_S, ctx);
//----------------------------------------------------
    return res;
}

void log_GMI(GArray* F, GArray* P, const Polynom h, int t, PolynomRing ctx){
    GArray* _P;
    GArray* rem_items;
    Polynom* ph;
    SPair* pspair;
    Polynom f, g;
    fq_nmod_mpoly_t lcm, div, lt_h, L, lt_f, gcd;
    int i, j, flag1, flag2, flag3;
//----------------------------------------------------
    _P = g_array_new(FALSE, FALSE, sizeof(SPair));
    rem_items = g_array_new(FALSE, FALSE, sizeof(ulong));
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_init(div, ctx);
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(L, ctx);
    fq_nmod_mpoly_init(lt_f, ctx);
    fq_nmod_mpoly_init(gcd, ctx);

    LT(lt_h, h, ctx);
//----------------------------------------------------
    for(i = 0; i < t; i++){
        SPair sp = {NULL, i, t};
        g_array_append_val(_P, sp);
    }

    i = 0;
    while(i < P->len){
        SPair sp = g_array_index(P, SPair, i);
        f = g_array_index(F, Polynom, sp.first);
        g = g_array_index(F, Polynom, sp.second);
        LCM(L, f, g, ctx);

        flag1 = fq_nmod_mpoly_divides(div, L, lt_h, ctx);

        LCM(lcm, h, f, ctx);
        flag2 = fq_nmod_mpoly_equal(lcm, L, ctx);

        LCM(lcm, h, g, ctx);
        flag3 = fq_nmod_mpoly_equal(lcm, L, ctx);

        if((flag1 == 1) && (flag2 == 0) && (flag3 == 0)){
            printf("rem pair(");
            fq_nmod_mpoly_print_pretty(sp.poly, NULL, ctx);
            printf(", %ld, %ld) from P\n", sp.first, sp.second);
            free_SPair(&sp, ctx);
            g_array_remove_index(P, i);
            i--;
        }
        i++;
    }

    printf("_P:\n");
    for(i = 0; i < _P->len; i++){
        printf("(NULL, %ld, %ld)\n", g_array_index(_P, SPair, i).first, g_array_index(_P, SPair, i).second);
    }

    i = 0;
    while(i < _P->len){
        printf("i: %d, _P->len: %d\n", i, _P->len);
        
        j = 0;
        while(j < _P->len){
            printf("j: %d\n", j);
            if (i != j){
                f = g_array_index(F, Polynom, g_array_index(_P, SPair, i).first);
                g = g_array_index(F, Polynom, g_array_index(_P, SPair, j).first);

                LCM(lcm, f, h, ctx);
                LCM(L, g, h, ctx);
                if (fq_nmod_mpoly_divides(div, L, lcm, ctx) == 1){
                    printf("rem pair: (%ld, %ld) from _P\n", g_array_index(_P, SPair, j).first, g_array_index(_P, SPair, j).second);
                    g_array_remove_index(_P, j);

                    if (j < i)
                        i--;
                }
            }
            j++;
        }
        i++;
    }

    printf("_P after:\n");
    for(i = 0; i < _P->len; i++){
        printf("(NULL, %ld, %ld)\n", g_array_index(_P, SPair, i).first, g_array_index(_P, SPair, i).second);
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, SPair, i).first);
        LT(lt_f, f, ctx);

        fq_nmod_mpoly_gcd(gcd, lt_f, lt_h, ctx);
        if (fq_nmod_mpoly_is_one(gcd, ctx) == 1){
            g_array_remove_index(_P, i);
        }
        else
            i++;
    }

    printf("_P after after:\n");
    for(i = 0; i < _P->len; i++){
        printf("(NULL, %ld, %ld)\n", g_array_index(_P, SPair, i).first, g_array_index(_P, SPair, i).second);
    }

    if (t == F->len){
        Polynom new_poly = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(new_poly, ctx);
        fq_nmod_mpoly_set(new_poly, h, ctx);

        g_array_append_val(F, new_poly);
    }

    pspair = (SPair*)_P->data;
    for(i = 0; i < _P->len; i++){
        SPair sp = {};
        f = g_array_index(F, Polynom, pspair->first);
        S(L, f, h, ctx); 

        if(fq_nmod_mpoly_is_zero(L, ctx) == 0){
            init_SPair(&sp, ctx);
            set_SPair(&sp, L, pspair->first, t, ctx);
            g_array_append_val(P, sp);
        }
        // S(sp.poly, f, h, ctx);
        // sp.first = pspair->first;
        // // sp.second = pspair->second;
        // sp.second = t;

        pspair++;
    }

    printf("P res:\n");
    for(i = 0; i < P->len; i++){
        printf("(");
        fq_nmod_mpoly_print_pretty(g_array_index(P, SPair, i).poly, NULL, ctx);
        printf(", %ld, %ld)\n", g_array_index(P, SPair, i).first, g_array_index(P, SPair, i).second);
    }

//----------------------------------------------------
    g_array_free(_P, TRUE);
    g_array_free(rem_items, TRUE);
    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(div, ctx);
    fq_nmod_mpoly_clear(lt_h, ctx);
    fq_nmod_mpoly_clear(L, ctx);
    fq_nmod_mpoly_clear(lt_f, ctx);
    fq_nmod_mpoly_clear(gcd, ctx);
}

void GMI(GArray* F, GArray* P, const Polynom h, int t, PolynomRing ctx){
    GArray* _P;
    GArray* rem_items;
    Polynom* ph;
    SPair* pspair;
    Polynom f, g;
    fq_nmod_mpoly_t lcm, div, lt_h, L, lt_f, gcd;
    int i, j, flag1, flag2, flag3;
//----------------------------------------------------
    _P = g_array_new(FALSE, FALSE, sizeof(SPair));
    rem_items = g_array_new(FALSE, FALSE, sizeof(ulong));
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_init(div, ctx);
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(L, ctx);
    fq_nmod_mpoly_init(lt_f, ctx);
    fq_nmod_mpoly_init(gcd, ctx);

    LT(lt_h, h, ctx);
//----------------------------------------------------
    for(i = 0; i < t; i++){
        SPair sp = {NULL, i, t};
        g_array_append_val(_P, sp);
    }

    i = 0;
    while(i < P->len){
        SPair sp = g_array_index(P, SPair, i);
        f = g_array_index(F, Polynom, sp.first);
        g = g_array_index(F, Polynom, sp.second);
        LCM(L, f, g, ctx);

        // flag1 = fq_nmod_mpoly_divides(div, L, lt_h, ctx);

        // LCM(lcm, h, f, ctx);
        // flag2 = fq_nmod_mpoly_equal(lcm, L, ctx);

        // LCM(lcm, h, g, ctx);
        // flag3 = fq_nmod_mpoly_equal(lcm, L, ctx);

        // if((flag1 == 1) && (flag2 == 0) && (flag3 == 0)){
        //     free_SPair(&sp, ctx);
        //     g_array_remove_index(P, i);
        //     i--;
        // }

        if (fq_nmod_mpoly_divides(div, L, lt_h, ctx) == 1){
            LCM(lcm, h, f, ctx);
            if (fq_nmod_mpoly_equal(lcm, L, ctx) == 0){
                LCM(lcm, h, g, ctx);
                if (fq_nmod_mpoly_equal(lcm, L, ctx) == 0){
                    free_SPair(&sp, ctx);
                    g_array_remove_index(P, i);
                    i--;
                }
            }
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, SPair, i).first);
        LCM(lcm, f, h, ctx);
        j = 0;
        while(j < _P->len){
            if (i != j){
                
                g = g_array_index(F, Polynom, g_array_index(_P, SPair, j).first);
                LCM(L, g, h, ctx);
                if (fq_nmod_mpoly_divides(div, L, lcm, ctx) == 1){
                    g_array_remove_index(_P, j);

                    if (j < i)
                        i--;

                    // continue ?
                }
            }
            j++;
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, SPair, i).first);
        LT(lt_f, f, ctx);

        fq_nmod_mpoly_gcd(gcd, lt_f, lt_h, ctx);
        if (fq_nmod_mpoly_is_one(gcd, ctx) == 1){
            g_array_remove_index(_P, i);
        }
        else
            i++;
    }

    if (t == F->len){
        Polynom new_poly = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(new_poly, ctx);
        fq_nmod_mpoly_set(new_poly, h, ctx);

        g_array_append_val(F, new_poly);
    }

    pspair = (SPair*)_P->data;
    for(i = 0; i < _P->len; i++){
        SPair sp = {};
        f = g_array_index(F, Polynom, pspair->first);
        S(L, f, h, ctx); 

        if(fq_nmod_mpoly_is_zero(L, ctx) == 0){
            init_SPair(&sp, ctx);
            set_SPair(&sp, L, pspair->first, t, ctx);
            g_array_append_val(P, sp);
        }

        pspair++;
    }

//----------------------------------------------------
    g_array_free(_P, TRUE);
    g_array_free(rem_items, TRUE);
    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(div, ctx);
    fq_nmod_mpoly_clear(lt_h, ctx);
    fq_nmod_mpoly_clear(L, ctx);
    fq_nmod_mpoly_clear(lt_f, ctx);
    fq_nmod_mpoly_clear(gcd, ctx);
}

Buchberger_result log_buchberger_v2(const Basis basis, ulong t, const PolynomRing ctx){
    GArray *F, *P;
    Polynom* hp;
    Polynom S_poly;
    Basis Q, G;
    fq_nmod_mpoly_t reminder;
    SPair sp;
    int i;
//----------------------------------------------------
    F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(F, f);
    }

    P = g_array_new(FALSE, FALSE, sizeof(SPair));
    fq_nmod_mpoly_init(reminder, ctx);
//----------------------------------------------------
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        log_GMI(F, P, *hp, i, ctx);
        hp++;
    }

    while(P->len > 0){
        if (F->len % 100 == 0)
            printf("F->len: %d\n", F->len);

        i = log_find_min(P, ctx);
        sp = g_array_index(P, SPair, i);

        printf("curr P:\n");
        for(int k = 0; k<P->len; k++){
            SPair sp2 = g_array_index(P, SPair, k);
            printf("(");
            fq_nmod_mpoly_print_pretty(sp2.poly, NULL, ctx);
            printf(", %ld, %ld)\n", sp2.first, sp2.second);
        }


        print_poly("min:", sp.poly, NULL, ctx);
        printf("\n");

        Q = init_empty_basis(F->len, ctx);
        fq_nmod_mpoly_divrem_ideal(Q, reminder, sp.poly, (Polynom*)F->data, F->len, ctx);
        free_basis(Q, F->len, ctx);
        print_poly("S_min mod F:", reminder, NULL, ctx);

        free_SPair(&sp, ctx);
        g_array_remove_index(P, i);
        
        if (fq_nmod_mpoly_is_zero(reminder, ctx) == 0){
            log_GMI(F, P, reminder, F->len, ctx);
        }

        printf("curr basis(%d):\n", F->len);
        for(i = 0; i<F->len; i++){
            fq_nmod_mpoly_print_pretty(g_array_index(F, Polynom, i), NULL, ctx);
            printf("\n");
        }
    }
//----------------------------------------------------
    Basis res = from_garray(F);
    ulong len = F->len;
    Buchberger_result resres = {res, len};
    fq_nmod_mpoly_clear(reminder, ctx);
    g_array_free(F, TRUE);
    g_array_free(P, TRUE);

    return resres;
}

Buchberger_result buchberger_v2(const Basis basis, ulong t, const PolynomRing ctx){
    GArray *F, *P;
    Polynom* hp;
    Polynom S_poly;
    // Basis Q, G;
    GArray* Q;
    fq_nmod_mpoly_t reminder;
    SPair sp;
    int i;
//----------------------------------------------------
    Q = g_array_new(FALSE, FALSE, sizeof(Polynom));
    F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(F, f);

        Polynom g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(g, ctx);
        g_array_append_val(Q, g);
    }

    P = g_array_new(FALSE, FALSE, sizeof(SPair));
    fq_nmod_mpoly_init(reminder, ctx);

    
    
//----------------------------------------------------
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        GMI(F, P, *hp, i, ctx);
        hp++;
    }

    while(P->len > 0){
        i = find_min(P, ctx);
        sp = g_array_index(P, SPair, i);

        // Q = init_empty_basis(F->len, ctx);
        fq_nmod_mpoly_divrem_ideal((Polynom*)Q->data, reminder, sp.poly, (Polynom*)F->data, F->len, ctx);
        // free_basis(Q, F->len, ctx);

        free_SPair(&sp, ctx);
        g_array_remove_index(P, i);
        
        if (fq_nmod_mpoly_is_zero(reminder, ctx) == 0){
            GMI(F, P, reminder, F->len, ctx);

            Polynom g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
            fq_nmod_mpoly_init(g, ctx);
            g_array_append_val(Q, g);
        }
    }
//----------------------------------------------------
    Basis res = from_garray(F);
    ulong len = F->len;
    Buchberger_result resres = {res, len};
    fq_nmod_mpoly_clear(reminder, ctx);
    g_array_free(F, TRUE);
    g_array_free(P, TRUE);

    for(i = 0; i < Q->len; i++){
        Polynom g = g_array_index(Q, Polynom, i);
        fq_nmod_mpoly_clear(g, ctx);
    }

    g_array_free(Q, TRUE);

    return resres;
}

int find_min_v2(GArray* F, GArray* P, PolynomRing ctx){
    fq_nmod_mpoly_t lcm_min, lcm_cur;
    Polynom f, g;
    int res;
    Pair* ppair;
    int i;
//----------------------------------------------------
    fq_nmod_mpoly_init(lcm_min, ctx);
    fq_nmod_mpoly_init(lcm_cur, ctx);
//----------------------------------------------------
    if (P->len > 0){
        ppair = (Pair*)P->data;
        f = g_array_index(F, Polynom, ppair->first);
        g = g_array_index(F, Polynom, ppair->second);

        LCM(lcm_min, f, g, ctx);
        res = 0;
        ppair++;

        for(i = 1; i < P->len; i++){
            f = g_array_index(F, Polynom, ppair->first);
            g = g_array_index(F, Polynom, ppair->second);
            LCM(lcm_cur, f, g, ctx);
            if (fq_nmod_mpoly_cmp(lcm_min, lcm_cur, ctx) == 1){
                fq_nmod_mpoly_set(lcm_min, lcm_cur, ctx);
                res = i;
            }
            ppair++;
        }
    }
//----------------------------------------------------
    fq_nmod_mpoly_clear(lcm_min, ctx);
    fq_nmod_mpoly_clear(lcm_cur, ctx);
//----------------------------------------------------
    return res;
}


void GMI_v2(GArray* F, GArray* P, const Polynom h, int t, PolynomRing ctx){
    GArray* _P;
    GArray* rem_items;
    Polynom* ph;
    Pair* pspair;
    Polynom f, g;
    fq_nmod_mpoly_t lcm, div, lt_h, L, lt_f, gcd;
    int i, j, flag1, flag2, flag3;
//----------------------------------------------------
    _P = g_array_new(FALSE, FALSE, sizeof(Pair));
    rem_items = g_array_new(FALSE, FALSE, sizeof(ulong));
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_init(div, ctx);
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(L, ctx);
    fq_nmod_mpoly_init(lt_f, ctx);
    fq_nmod_mpoly_init(gcd, ctx);

    LT(lt_h, h, ctx);
//----------------------------------------------------
    for(i = 0; i < t; i++){
        Pair sp = {i, t};
        g_array_append_val(_P, sp);
    }

    i = 0;
    while(i < P->len){
        Pair sp = g_array_index(P, Pair, i);
        f = g_array_index(F, Polynom, sp.first);
        g = g_array_index(F, Polynom, sp.second);
        LCM(L, f, g, ctx);

        if (fq_nmod_mpoly_divides(div, L, lt_h, ctx) == 1){
            LCM(lcm, h, f, ctx);
            if (fq_nmod_mpoly_equal(lcm, L, ctx) == 0){
                LCM(lcm, h, g, ctx);
                if (fq_nmod_mpoly_equal(lcm, L, ctx) == 0){
                    g_array_remove_index(P, i);
                    i--;
                }
            }
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, Pair, i).first);
        LCM(lcm, f, h, ctx);
        j = 0;
        while(j < _P->len){
            if (i != j){
                
                g = g_array_index(F, Polynom, g_array_index(_P, Pair, j).first);
                LCM(L, g, h, ctx);
                if (fq_nmod_mpoly_divides(div, L, lcm, ctx) == 1){
                    g_array_remove_index(_P, j);

                    if (j < i)
                        i--;

                    continue;
                }
            }
            j++;
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, Pair, i).first);
        LT(lt_f, f, ctx);

        fq_nmod_mpoly_gcd(gcd, lt_f, lt_h, ctx);
        if (fq_nmod_mpoly_is_one(gcd, ctx) == 1){
            g_array_remove_index(_P, i);
        }
        else
            i++;
    }

    if (t == F->len){
        Polynom new_poly = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(new_poly, ctx);
        fq_nmod_mpoly_set(new_poly, h, ctx);

        g_array_append_val(F, new_poly);
    }

    pspair = (Pair*)_P->data;
    for(i = 0; i < _P->len; i++){
        Pair sp = {pspair->first, t};
        g_array_append_val(P, sp);
        pspair++;
    }

//----------------------------------------------------
    g_array_free(_P, TRUE);
    g_array_free(rem_items, TRUE);
    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(div, ctx);
    fq_nmod_mpoly_clear(lt_h, ctx);
    fq_nmod_mpoly_clear(L, ctx);
    fq_nmod_mpoly_clear(lt_f, ctx);
    fq_nmod_mpoly_clear(gcd, ctx);
}

Buchberger_result buchberger_v2_1(const Basis basis, ulong t, const PolynomRing ctx){
    GArray *F, *P;
    Polynom* hp;
    Polynom a, b;
    GArray* Q;
    fq_nmod_mpoly_t reminder, S_poly;
    Pair sp;
    int i;
//----------------------------------------------------
    Q = g_array_new(FALSE, FALSE, sizeof(Polynom));
    F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(F, f);

        Polynom g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(g, ctx);
        g_array_append_val(Q, g);
    }

    P = g_array_new(FALSE, FALSE, sizeof(Pair));
    fq_nmod_mpoly_init(reminder, ctx);
    fq_nmod_mpoly_init(S_poly, ctx);

//----------------------------------------------------
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        GMI_v2(F, P, *hp, i, ctx);
        hp++;
    }

    while(P->len > 0){
        i = find_min_v2(F, P, ctx);
        sp = g_array_index(P, Pair, i);
        a = g_array_index(F, Polynom, sp.first);
        b = g_array_index(F, Polynom, sp.second);
        S(S_poly, a, b, ctx);

        // Q = init_empty_basis(F->len, ctx);
        fq_nmod_mpoly_divrem_ideal((Polynom*)Q->data, reminder, S_poly, (Polynom*)F->data, F->len, ctx);
        // free_basis(Q, F->len, ctx);

        g_array_remove_index(P, i);
        
        if (fq_nmod_mpoly_is_zero(reminder, ctx) == 0){
            GMI_v2(F, P, reminder, F->len, ctx);

            Polynom g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
            fq_nmod_mpoly_init(g, ctx);
            g_array_append_val(Q, g);
        }
    }
//----------------------------------------------------
    Basis res = from_garray(F);
    ulong len = F->len;
    Buchberger_result resres = {res, len};
    fq_nmod_mpoly_clear(reminder, ctx);
    g_array_free(F, FALSE);
    g_array_free(P, TRUE);

    for(i = 0; i < Q->len; i++){
        Polynom g = g_array_index(Q, Polynom, i);
        fq_nmod_mpoly_clear(g, ctx);
    }

    g_array_free(Q, TRUE);
    fq_nmod_mpoly_clear(S_poly, ctx);

    return resres;
}

int find_max_from_indexes(GArray* F, GArray* P, GArray* indx, PolynomRing ctx){
    fq_nmod_mpoly_t lcm_max, lcm_cur;
    Polynom f, g;
    int res;
    Pair pair;
    int i, j;

//----------------------------------------------------
    fq_nmod_mpoly_init(lcm_max, ctx);
    fq_nmod_mpoly_init(lcm_cur, ctx);
//----------------------------------------------------
    // printf("------------------------find_max-------------------------\n");
    if (P->len > 0 && indx->len > 0){
        pair = g_array_index(P, Pair, g_array_index(indx, ulong, 0));
        f = g_array_index(F, Polynom, pair.first);
        g = g_array_index(F, Polynom, pair.second);

        LCM(lcm_max, f, g, ctx);
        res = 0;

        // print_poly("lcm_max", lcm_max, NULL, ctx);
        

        for(i = 1; i < indx->len; i++){
            pair = g_array_index(P, Pair, g_array_index(indx, ulong, i));
            f = g_array_index(F, Polynom, pair.first);
            g = g_array_index(F, Polynom, pair.second);
            LCM(lcm_cur, f, g, ctx);
            // print_poly("lcm_cur", lcm_cur, NULL, ctx);
            // printf("lcm_cur > lcm_max: %d\n", fq_nmod_mpoly_cmp(lcm_cur, lcm_max, ctx));
            if (fq_nmod_mpoly_cmp(lcm_cur, lcm_max, ctx) == 1){
                fq_nmod_mpoly_set(lcm_max, lcm_cur, ctx);
                res = i;
            }
        }
    }
    // printf("--------------------------------------------------------\n");
//----------------------------------------------------
    fq_nmod_mpoly_clear(lcm_max, ctx);
    fq_nmod_mpoly_clear(lcm_cur, ctx);
//----------------------------------------------------
    return res;
}

void find_min_v3(GArray* F, GArray* P, ulong* res, ulong count, PolynomRing ctx){
    fq_nmod_mpoly_t lcm_min, lcm_cur, lcm_min_max;
    Polynom f, g;
    Pair pair;
    ulong i, j, min_max;

    GArray* min_elems = g_array_new(FALSE, FALSE, sizeof(ulong));
//----------------------------------------------------
    fq_nmod_mpoly_init(lcm_min, ctx);
    fq_nmod_mpoly_init(lcm_cur, ctx);
    fq_nmod_mpoly_init(lcm_min_max, ctx);
    
    for(i = 0; i < count; i++){
        res[i] = -1;
    }

    if(count >= P->len){
        for(i = 0; i < P->len; i++){
            res[i] = i;
        }
        return;
    }
//----------------------------------------------------
    if (P->len > 0){
        for(i = 0; i < count; i++)
            g_array_append_val(min_elems, i);

        min_max = find_max_from_indexes(F, P, min_elems, ctx);
        // printf("min_max=%ld\n", min_max);
        pair = g_array_index(P, Pair, g_array_index(min_elems, ulong , min_max));
        f = g_array_index(F, Polynom, pair.first);
        g = g_array_index(F, Polynom, pair.second);
        LCM(lcm_min_max, f, g, ctx);
        

        for(i = count; i < P->len; i++){
            // printf("min_elems:");
            // print_ulong_garray(min_elems);

            // for(int k = 0; k < min_elems->len; k++){
            //     pair = g_array_index(P, Pair, g_array_index(min_elems, ulong, k));
            //     f = g_array_index(F, Polynom, pair.first);
            //     g = g_array_index(F, Polynom, pair.second);
            //     LCM(lcm_cur, f, g, ctx);
            //     print_poly("", lcm_cur, NULL, ctx);
            // }

            pair = g_array_index(P, Pair, i);
            f = g_array_index(F, Polynom, pair.first);
            g = g_array_index(F, Polynom, pair.second);
            LCM(lcm_cur, f, g, ctx);
            // printf("(%ld, %ld)\n", pair.first, pair.second);
            // print_poly("lcm_min_max:", lcm_min_max, NULL, ctx);
            // print_poly("lcm_cur:", lcm_cur, NULL, ctx);
            // printf("lcm_min_max>lcm_cur: %d\n", fq_nmod_mpoly_cmp(lcm_min_max, lcm_cur, ctx));
            if (fq_nmod_mpoly_cmp(lcm_min_max, lcm_cur, ctx) == 1){
                ((ulong*)min_elems->data)[min_max] = i;
                min_max = find_max_from_indexes(F, P, min_elems, ctx);
                // printf("min_max=%ld\n", min_max);
                
                pair = g_array_index(P, Pair, g_array_index(min_elems, ulong , min_max));
                f = g_array_index(F, Polynom, pair.first);
                g = g_array_index(F, Polynom, pair.second);
                LCM(lcm_min_max, f, g, ctx);
            }
        }
    }

    for(i = 0; i < count; i++)
        res[i] = g_array_index(min_elems, ulong, i);
//----------------------------------------------------
    fq_nmod_mpoly_clear(lcm_min, ctx);
    fq_nmod_mpoly_clear(lcm_cur, ctx);
    fq_nmod_mpoly_clear(lcm_min_max, ctx);
    g_array_free(min_elems, TRUE);
}

int check_cond(thread_buchberger_data_t* data){
    return (*(data->finished_threads) >= NO_OF_IRRED);
}

// is monom f divide monom g
// int is_monom_divide(Polynom f, Polynom g, PolynomRing ctx){
//     ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
//     ulong exp_f[nvars];
//     ulong exp_g[nvars];
//     ulong exp_monom[nvars];


//     fq_nmod_mpoly_get_term_exp_ui(exp_f, f, 0, ctx);
//     fq_nmod_mpoly_get_term_exp_ui(exp_g, g, 0, ctx);

//     for(int i = 0; i < nvars; i++){
//         if (exp_f[i] > exp_g[i])
//             return 0;
//     }

//     return 1;
// }

void* log_buchberger_task(void* params){
    thread_buchberger_data_t* data = (thread_buchberger_data_t*)params;
    fq_nmod_mpoly_t reminder, lt_f, lt_p, div, mul, cpy;
    ulong i, isdiv;

    if (check_cond(data) == 1) my_exit(data);
    if (fq_nmod_mpoly_is_zero(data->S_poly, data->ctx) == 1) pthread_exit(NULL);

    printf("--------------------div--------------------\n");
    // print_poly("S_poly", data->S_poly, NULL, data->ctx);

//-------------------------------------------------------------------------------
    fq_nmod_mpoly_init(reminder, data->ctx);
    fq_nmod_mpoly_zero(reminder, data->ctx);
    fq_nmod_mpoly_init(lt_f, data->ctx);
    fq_nmod_mpoly_init(lt_p, data->ctx);
    fq_nmod_mpoly_init(div, data->ctx);
    fq_nmod_mpoly_init(mul, data->ctx);
    fq_nmod_mpoly_init(cpy, data->ctx);
//-------------------------------------------------------------------------------
    while (fq_nmod_mpoly_is_zero(data->S_poly, data->ctx) == 0){
        if (check_cond(data) == 1) {
            printf("pthread_exit cond is true");
            my_exit(data);
        }
        print_poly("S_poly:", data->S_poly, NULL, data->ctx);
        i = 0;
        isdiv = 0;
        LT(lt_p, data->S_poly, data->ctx);
        print_poly("LT(p):", lt_p, NULL, data->ctx);
        while((i < data->npoli) && (isdiv==0)){
            if (check_cond(data) == 1) {
                printf("pthread_exit cond is true");
                my_exit(data);
            }
            LT(lt_f, data->F[i], data->ctx);
            print_poly("LT(f):", lt_f, NULL, data->ctx);
            printf("npoli=%ld, i=%ld, isdiv=%ld, %d\n", data->npoli, i, isdiv, fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx));
            if (fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx) == 1){
                fq_nmod_mpoly_mul(mul, div, data->F[i], data->ctx);
                fq_nmod_mpoly_set(cpy, data->S_poly, data->ctx);
                fq_nmod_mpoly_sub(data->S_poly, cpy, mul, data->ctx);
                
                isdiv = 1;
            } else{
                printf("false\n");
                i++;
                printf("npoli=%ld, i=%ld, isdiv=%ld, %d\n", data->npoli, i, isdiv, fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx));
            }
        }
        if (isdiv == 0){
            // print_poly("lt_p:", lt_p, NULL, data->ctx);
            fq_nmod_mpoly_set(cpy, reminder, data->ctx);
            fq_nmod_mpoly_add(reminder, cpy, lt_p, data->ctx);
            
            fq_nmod_mpoly_set(cpy, data->S_poly, data->ctx);
            fq_nmod_mpoly_sub(data->S_poly, cpy, lt_p, data->ctx);
        }
    }

    fq_nmod_mpoly_set(data->S_poly, reminder, data->ctx);
//-------------------------------------------------------------------------------
    fq_nmod_mpoly_clear(reminder, data->ctx);
    fq_nmod_mpoly_clear(lt_f, data->ctx);
    fq_nmod_mpoly_clear(lt_p, data->ctx);
    fq_nmod_mpoly_clear(div, data->ctx);
    fq_nmod_mpoly_clear(mul, data->ctx);
    fq_nmod_mpoly_clear(cpy, data->ctx);

    data->completed = 1;
}

void my_exit(thread_buchberger_data_t* data){
    fq_nmod_mpoly_zero(data->S_poly, data->ctx);
    data->completed = 0;
    pthread_exit(NULL);
}

void* buchberger_task(void* params){
    thread_buchberger_data_t* data = (thread_buchberger_data_t*)params;
    fq_nmod_mpoly_t reminder, lt_f, lt_p, div, mul, cpy;
    ulong i, isdiv;

    if (check_cond(data) == 1) my_exit(data);
    if (fq_nmod_mpoly_is_zero(data->S_poly, data->ctx) == 1) pthread_exit(NULL);

    // printf("--------------------div--------------------\n");
    // print_poly("S_poly", data->S_poly, NULL, data->ctx);

//-------------------------------------------------------------------------------
    fq_nmod_mpoly_init(reminder, data->ctx);
    fq_nmod_mpoly_zero(reminder, data->ctx);
    fq_nmod_mpoly_init(lt_f, data->ctx);
    fq_nmod_mpoly_init(lt_p, data->ctx);
    fq_nmod_mpoly_init(div, data->ctx);
    fq_nmod_mpoly_init(mul, data->ctx);
    fq_nmod_mpoly_init(cpy, data->ctx);
//-------------------------------------------------------------------------------
    while (fq_nmod_mpoly_is_zero(data->S_poly, data->ctx) == 0){
        // print_poly("S_poly:", data->S_poly, NULL, data->ctx);
        if (check_cond(data) == 1) my_exit(data);
        i = 0;
        isdiv = 0;
        LT(lt_p, data->S_poly, data->ctx);
        // print_poly("LT(p):", lt_p, NULL, data->ctx);
        while((i < data->npoli) && (isdiv==0)){
            if (check_cond(data) == 1) my_exit(data);
            LT(lt_f, data->F[i], data->ctx);
            // print_poly("LT(f):", lt_f, NULL, data->ctx);
            // printf("npoli=%d, i=%d, isdiv=%d, %d\n", data->npoli, i, isdiv, fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx));
            if (fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx) == 1){
                fq_nmod_mpoly_mul(mul, div, data->F[i], data->ctx);
                fq_nmod_mpoly_set(cpy, data->S_poly, data->ctx);
                fq_nmod_mpoly_sub(data->S_poly, cpy, mul, data->ctx);
                
                isdiv = 1;
            } else{
                // printf("false\n");
                i++;
                // printf("npoli=%d, i=%d, isdiv=%d, %d\n", data->npoli, i, isdiv, fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx));
            }
        }
        if (isdiv == 0){
            // print_poly("lt_p:", lt_p, NULL, data->ctx);
            fq_nmod_mpoly_set(cpy, reminder, data->ctx);
            fq_nmod_mpoly_add(reminder, cpy, lt_p, data->ctx);
            
            fq_nmod_mpoly_set(cpy, data->S_poly, data->ctx);
            fq_nmod_mpoly_sub(data->S_poly, cpy, lt_p, data->ctx);
        }
    }

    fq_nmod_mpoly_set(data->S_poly, reminder, data->ctx);
//-------------------------------------------------------------------------------
    fq_nmod_mpoly_clear(reminder, data->ctx);
    fq_nmod_mpoly_clear(lt_f, data->ctx);
    fq_nmod_mpoly_clear(lt_p, data->ctx);
    fq_nmod_mpoly_clear(div, data->ctx);
    fq_nmod_mpoly_clear(mul, data->ctx);
    fq_nmod_mpoly_clear(cpy, data->ctx);

    data->completed = 1;
}

int find_min_reduction(Basis reminders, ulong n, PolynomRing ctx){
    Polynom min;
    int i = 0;
    int res = -1;
    while(i < n){
        if (fq_nmod_mpoly_is_zero(reminders[i], ctx) == 0){
            min = reminders[i];
            res = i;
            i++;
            break;
        }
        i++;
    }

    for (i; i < n; i++){
        if (fq_nmod_mpoly_is_zero(reminders[i], ctx) == 0 & fq_nmod_mpoly_cmp(min, reminders[i], ctx) == 1){
            min = reminders[i];
            res = i;
        }
    }

    return res;
}

int find_min_reduction_v2(thread_buchberger_data_v2_t* datas, ulong n, PolynomRing ctx){
    Polynom min;
    int i = 0;
    int res = -1;
    while(i < n){
        if (datas[i].executed == 1 && fq_nmod_mpoly_is_zero(datas[i].S_poly, ctx) == 0){
            min = datas[i].S_poly;
            res = i;
            i++;
            break;
        }
        i++;
    }

    // printf("res=%d\n", res);

    for (i; i < n; i++){
        if ( datas[i].executed == 1 && fq_nmod_mpoly_is_zero(datas[i].S_poly, ctx) == 0 && fq_nmod_mpoly_cmp(min, datas[i].S_poly, ctx) == 1){
            min = datas[i].S_poly;
            res = i;
        }
    }

    return res;
}

Buchberger_result log_threaded_buchberger(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx){
    GArray *F, *P;
    Polynom* hp;
    Basis S_polys;
    thread_buchberger_data_t* data;
    pthread_t* threads;
    pthread_attr_t* attrs;
    Pair sp;
    Polynom f, g;
    ulong i, j, k;
    ulong finished_threads;
    ulong* selected_pairs;
    Pair* buffer = calloc(threads_count, sizeof(Pair));

//----------------------------------------------------
    F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(F, f);
    }

    P = g_array_new(FALSE, FALSE, sizeof(Pair));

    data = calloc(threads_count, sizeof(thread_buchberger_data_t));
    threads = calloc(threads_count, sizeof(pthread_t));
    attrs = calloc(threads_count, sizeof(pthread_attr_t));

    S_polys = init_empty_basis(threads_count, ctx);

    for(i = 0; i < threads_count; i++){
        pthread_attr_init(&attrs[i]);
        Polynom S_poly, reminder;

        data[i].finished_threads = &finished_threads;
        data[i].S_poly = S_polys[i];
        data[i].ctx = ctx;
    }

    selected_pairs = flint_calloc(threads_count, sizeof(ulong));

//----------------------------------------------------
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        GMI_v2(F, P, *hp, i, ctx);
        hp++;
    }

    g_array_sort(P, cmpPair);

    printf("curr P(%d):\n", P->len);
    for(i = 0; i<P->len; i++){
        Pair sp2 = g_array_index(P, Pair, i);
        printf("(%ld, %ld)\n", sp2.first, sp2.second);
    }

    while(P->len > 0){
        j = 0;

        fq_nmod_mpoly_t lcm;
        fq_nmod_mpoly_init(lcm, ctx);

        find_min_v3(F, P, selected_pairs, threads_count, ctx);
        printf("Selected pairs:\n");
        j = 0;
        for(i = 0; i < threads_count; i++){
            k = selected_pairs[i];
            if (k != -1){
                sp = g_array_index(P, Pair, k);
                f = g_array_index(F, Polynom, sp.first);
                g = g_array_index(F, Polynom, sp.second);
                
                S(S_polys[i], f, g, ctx);
                // LCM(lcm, f, g, ctx);
                // print_poly("lcm:", lcm, NULL, ctx);
                printf("(%ld, %ld)\n", sp.first, sp.second);
                j++;
            }
            else{
                printf("None\n");
            }
        }
            
        

        // printf("Selected pairs:\n");
        // while(j < threads_count && P->len > 0){
        //     k = find_min_v2(F, P, ctx);
        //     sp = g_array_index(P, Pair, k);
        //     f = g_array_index(F, Polynom, sp.first);
        //     g = g_array_index(F, Polynom, sp.second);
        //     g_array_remove_index(P, k);
            
        //     // LCM(lcm, f, g, ctx);
        //     // print_poly("lcm:", lcm, NULL, ctx);

        //     printf("(%ld, %ld)\n", sp.first, sp.second);

        //     S(S_polys[j], f, g, ctx);
        //     j++;
        // }

        // fq_nmod_mpoly_clear(lcm, ctx);

        if (j < threads_count){
            fq_nmod_mpoly_zero(S_polys[j], ctx);
            j++;
        }

        printf("S_polys:\n");
        for(i = 0; i < threads_count; i++){
            printf("i=%ld: ", i);
            fq_nmod_mpoly_print_pretty(S_polys[i], NULL, ctx);
            printf("\n");
        }


        // j = 0;
        finished_threads = 0;
        for (i = 0; i < threads_count; i++){
            data[i].F = (Polynom*)F->data;
            data[i].npoli = F->len;
            data[i].completed = 0;
            pthread_create(&threads[i], &attrs[i], buchberger_task, &data[i]);
        }

        for (i = 0; i < threads_count; i++)
            pthread_join(threads[i], NULL);

        printf("reminders:\n");
        for(i = 0; i < threads_count; i++){
            printf("i=%ld: ", i);
            fq_nmod_mpoly_print_pretty(S_polys[i], NULL, ctx);
            printf("\n");
        }

        k = find_min_reduction(S_polys, threads_count, ctx);
        printf("Selected k=%ld\n", k);

        if (k != -1){
            g_array_remove_index(P, selected_pairs[k]);
            GMI_v2(F, P, S_polys[k], F->len, ctx);
            g_array_sort(P, cmpPair);
        } else if (P->len > 0){
            // for(i = 0; i < threads_count; i++){
            //     // if (selected_pairs[i] != -1){
            //     //     // Pair pair = g_array_index(P, Pair, selected_pairs[i]);
            //     //     // printf("Deleition pair %ld: (%ld, %ld)\n", selected_pairs[i], pair.first, pair.second);
            //     //     // printf("curr P(%d):\n", P->len);
            //     //     // for(int kk = 0; kk<P->len; kk++){
            //     //     //     Pair sp2 = g_array_index(P, Pair, kk);
            //     //     //     printf("(%ld, %ld)\n", sp2.first, sp2.second);
            //     //     // }

            //     //     // printf("index=%ld, P->len=%d\n", selected_pairs[i], P->len);
            //     //     // g_array_remove_index(P, selected_pairs[i]);

            //     //     Pair pair = g_array_index(P, Pair, selected_pairs[i]);
            //     //     int ind = 0;
            //     //     g_array_binary_search(P, &pair, cmpPair, &ind);
            //     //     g_array_remove_index(P, ind);
            //     // }
                
            //     // for(j = 0; j < threads_count; j++){
            //     //     Pair pair = g_array_index(P, Pair, selected_pairs[j]);
            //     //     // printf("(%ld, %ld)\n", pair.first, pair.second);
            //     //     printf("i=%ld, selectd_pairs[%ld]=%ld (%ld, %ld)\n", j, j, selected_pairs[j], pair.first, pair.second);
            //     // }
            //     // if (selected_pairs[i] != -1 && data[i].completed == 1){
            //     //     // g_array_remove_index(P, selected_pairs[i]);
                    
            //     //     printf("P(%d):\n", P->len);
            //     //     for(j = 0; j<P->len; j++){
            //     //         Pair sp2 = g_array_index(P, Pair, j);
            //     //         printf("(%ld, %ld)\n", sp2.first, sp2.second);
            //     //     }

            //     //     Pair pair = g_array_index(P, Pair, selected_pairs[i]);
            //     //     printf("remove pair (%ld, %ld)\n", pair.first, pair.second);
            //     //     int ind = -1;
                    
            //     //     g_array_binary_search(P, &pair, cmpPair, &ind);
            //     //     printf("P->len=%d, ind=%d\n", P->len, ind);
            //     //     g_array_remove_index(P, ind);
            //     // 
            // }

            ulong threads_completed = 0;
            for(j = 0; j < threads_count; j++)
                if (selected_pairs[j] != -1)
                    threads_completed++;

            k = 0;
            Pair pair;
            for(j = 0; j < threads_count; j++)
                if (selected_pairs[j] != -1){
                    pair = g_array_index(P, Pair, selected_pairs[j]);
                    buffer[k].first = pair.first;
                    buffer[k].second = pair.second;
                    k++;
                }

            int idx = 0;
            for(j = 0; j < threads_completed; j++){
                g_array_binary_search(P, &buffer[j], cmpPair, &idx);
                g_array_remove_index(P, idx);
            }

        }

        printf("curr P(%d):\n", P->len);
        for(i = 0; i<P->len; i++){
            Pair sp2 = g_array_index(P, Pair, i);
            printf("(%ld, %ld)\n", sp2.first, sp2.second);
        }

        printf("curr basis(%d):\n", F->len);
        // for(i = 0; i<F->len; i++){
        //     fq_nmod_mpoly_print_pretty(g_array_index(F, Polynom, i), NULL, ctx);
        //     printf("\n");
        // }

        // sleep(1);        
        // break;
    }
//----------------------------------------------------
    Basis res = from_garray(F);
    ulong len = F->len;
    Buchberger_result resres = {res, len};
    // fq_nmod_mpoly_clear(reminder, ctx);
    g_array_free(F, TRUE);
    g_array_free(P, TRUE);

    for(i = 0; i < threads_count; i++){
        pthread_attr_destroy(&attrs[i]);
    }
    free_basis(S_polys, threads_count, ctx);
    // fq_nmod_mpoly_clear(S_poly, ctx);
    free(buffer);

    return resres;
}

Buchberger_result threaded_buchberger(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx){
    GArray *F, *P;
    Polynom* hp;
    Basis S_polys;
    thread_buchberger_data_t* data;
    pthread_t* threads;
    pthread_attr_t* attrs;
    Pair sp;
    Polynom f, g;
    int i, j, k;
    ulong finished_threads;
    ulong* selected_pairs;
    Pair* buffer = calloc(threads_count, sizeof(Pair));

//----------------------------------------------------
    F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(F, f);
    }

    P = g_array_new(FALSE, FALSE, sizeof(Pair));

    data = calloc(threads_count, sizeof(thread_buchberger_data_t));
    threads = calloc(threads_count, sizeof(pthread_t));
    attrs = calloc(threads_count, sizeof(pthread_attr_t));

    S_polys = init_empty_basis(threads_count, ctx);

    for(i = 0; i < threads_count; i++){
        pthread_attr_init(&attrs[i]);
        Polynom S_poly, reminder;

        data[i].finished_threads = &finished_threads;
        data[i].S_poly = S_polys[i];
        data[i].ctx = ctx;
    }

    selected_pairs = flint_calloc(threads_count, sizeof(ulong));

//----------------------------------------------------
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        GMI_v2(F, P, *hp, i, ctx);
        hp++;
    }

    g_array_sort(P, cmpPair);

    while(P->len > 0){
        j = 0;

        fq_nmod_mpoly_t lcm;
        fq_nmod_mpoly_init(lcm, ctx);

        find_min_v3(F, P, selected_pairs, threads_count, ctx);

        j = 0;
        for(i = 0; i < threads_count; i++){
            k = selected_pairs[i];
            if (k != -1){
                sp = g_array_index(P, Pair, k);
                f = g_array_index(F, Polynom, sp.first);
                g = g_array_index(F, Polynom, sp.second);
                
                S(S_polys[i], f, g, ctx);
                j++;
            }
        }
            
        if (j < threads_count){
            // fq_nmod_mpoly_zero(S_polys[j], ctx);
            // j++;
            for(j; j < threads_count; j++)
                fq_nmod_mpoly_zero(S_polys[j], ctx);
        }


        // j = 0;
        finished_threads = 0;
        for (i = 0; i < threads_count; i++){
            
            if (fq_nmod_mpoly_is_zero(data[i].S_poly, ctx) == 1) continue;
            // print_poly("S_poly:", data[i].S_poly, NULL, ctx);
            
            data[i].F = (Polynom*)F->data;
            data[i].npoli = F->len;
            data[i].completed = 0;
            pthread_create(&threads[i], &attrs[i], buchberger_task, &data[i]);
        }

        for (i = 0; i < threads_count; i++){
            if (fq_nmod_mpoly_is_zero(data[i].S_poly, ctx) == 1) continue;
            pthread_join(threads[i], NULL);
        }


        k = find_min_reduction(S_polys, threads_count, ctx);

        if (k != -1){
            g_array_remove_index(P, selected_pairs[k]);
            GMI_v2(F, P, S_polys[k], F->len, ctx);
        } else if (P->len > 0){
            // g_array_sort(P, cmpPair);

            // ulong threads_completed = 0;
            // for(j = 0; j < threads_count; j++)
            //     if (selected_pairs[j] != -1)
            //         threads_completed++;

            // k = 0;
            // Pair pair;
            // for(j = 0; j < threads_count; j++)
            //     if (selected_pairs[j] != -1){
            //         pair = g_array_index(P, Pair, selected_pairs[j]);
            //         buffer[k].first = pair.first;
            //         buffer[k].second = pair.second;
            //         k++;
            //     }

            // int idx = 0;
            // for(j = 0; j < threads_completed; j++){
            //     g_array_binary_search(P, &buffer[j], cmpPair, &idx);
            //     g_array_remove_index(P, idx);
            // }

            // for(j = threads_count-1; j >= 0; j--){
            //     Pair pair = g_array_index(P, Pair, selected_pairs[j]);
            //     // printf("j=%d, %d\n", j, (j>=0));
            //     printf("i=%d, selectd_pairs[%d]=%ld P->len=%d (%ld, %ld)\n", j, j, selected_pairs[j], P->len, pair.first, pair.second);
            // }

            
            quick_sort(selected_pairs, 0, threads_count-1);

            // printf("after:\n");

            // for(j = threads_count-1; j >= 0; j--){
            //     Pair pair = g_array_index(P, Pair, selected_pairs[j]);
            //     // printf("j=%d, %d\n", j, (j>=0));
            //     printf("i=%d, selectd_pairs[%d]=%ld P->len=%d (%ld, %ld)\n", j, j, selected_pairs[j], P->len, pair.first, pair.second);
            // }

            // printf("hello\n");

            for(j = threads_count-1; j >= 0; j--)
                if(selected_pairs[j] != -1)
                    g_array_remove_index(P, selected_pairs[j]);

        }

    }
//----------------------------------------------------
    Basis res = from_garray(F);
    ulong len = F->len;
    Buchberger_result resres = {res, len};
    // fq_nmod_mpoly_clear(reminder, ctx);
    g_array_free(F, TRUE);
    g_array_free(P, TRUE);

    for(i = 0; i < threads_count; i++){
        pthread_attr_destroy(&attrs[i]);
    }
    free_basis(S_polys, threads_count, ctx);
    free(data);
    free(threads);
    free(attrs);
    // fq_nmod_mpoly_clear(S_poly, ctx);
    free(buffer);
    flint_free(selected_pairs);

    return resres;
}

void* buchberger_task_v2(void* params){
    thread_buchberger_data_v2_t* data = (thread_buchberger_data_v2_t*)params;
    fq_nmod_mpoly_t reminder, lt_f, lt_p, div, mul, cpy;
    ulong i, isdiv;
//-------------------------------------------------------------------------------
    fq_nmod_mpoly_init(reminder, data->ctx);
    fq_nmod_mpoly_zero(reminder, data->ctx);
    fq_nmod_mpoly_init(lt_f, data->ctx);
    fq_nmod_mpoly_init(lt_p, data->ctx);
    fq_nmod_mpoly_init(div, data->ctx);
    fq_nmod_mpoly_init(mul, data->ctx);
    fq_nmod_mpoly_init(cpy, data->ctx);
//-------------------------------------------------------------------------------
    while (1){
        loop1:
        // printf("%d %d\n", data->running, data->pause);
        if (data->running == 0) break;
        if (data->pause == 1) continue;

        // if (data->S_poly == NULL) continue;
        // print_poly("LOL", data->S_poly, NULL, data->ctx);
        while (fq_nmod_mpoly_is_zero(data->S_poly, data->ctx) == 0){
            // print_poly("S_poly:", data->S_poly, NULL, data->ctx);
            i = 0;
            isdiv = 0;
            LT(lt_p, data->S_poly, data->ctx);
            while((i < data->F->len) && (isdiv==0)){
                if (*(data->executed_threads) >= NO_OF_IRRED){
                    data->pause = 1;
                    goto loop1;
                }
                LT(lt_f, g_array_index(data->F, Polynom, i), data->ctx);
                if (fq_nmod_mpoly_divides(div, lt_p, lt_f, data->ctx) == 1){
                    fq_nmod_mpoly_mul(mul, div, g_array_index(data->F, Polynom, i), data->ctx);
                    fq_nmod_mpoly_set(cpy, data->S_poly, data->ctx);
                    fq_nmod_mpoly_sub(data->S_poly, cpy, mul, data->ctx);
                    isdiv = 1;
                } else{
                    i++;
                }
            }
            if (isdiv == 0){
                // print_poly("lt_p:", lt_p, NULL, data->ctx);
                fq_nmod_mpoly_set(cpy, reminder, data->ctx);
                fq_nmod_mpoly_add(reminder, cpy, lt_p, data->ctx);
                
                fq_nmod_mpoly_set(cpy, data->S_poly, data->ctx);
                fq_nmod_mpoly_sub(data->S_poly, cpy, lt_p, data->ctx);
            }
        }

        // print_poly("", reminder, NULL, data->ctx);
        // printf("YEY\n");
        fq_nmod_mpoly_set(data->S_poly, reminder, data->ctx);
        *(data->executed_threads)+=1;
        // printf("finished_threads=%d\n", *(data->finished_threads));
        data->executed = 1;
        data->pause = 1;
        fq_nmod_mpoly_zero(reminder, data->ctx);
        
    }
//-------------------------------------------------------------------------------
    fq_nmod_mpoly_clear(reminder, data->ctx);
    fq_nmod_mpoly_clear(lt_f, data->ctx);
    fq_nmod_mpoly_clear(lt_p, data->ctx);
    fq_nmod_mpoly_clear(div, data->ctx);
    fq_nmod_mpoly_clear(mul, data->ctx);
    fq_nmod_mpoly_clear(cpy, data->ctx);
}

void GMI_v3(GArray* F, GArray* P, const Polynom h, int t, thread_buchberger_data_v2_t* data, int data_len, PolynomRing ctx){
    GArray* _P;
    GArray* rem_items;
    Polynom* ph;
    Pair* pspair;
    Polynom f, g;
    fq_nmod_mpoly_t lcm, div, lt_h, L, lt_f, gcd;
    int i, j, flag1, flag2, flag3;
    int flag;

    if (h == NULL){
        i = 0;
        while(i < P->len){
            flag = 0;
            Pair sp = g_array_index(P, Pair, i);
            for(j = 0; j < data_len; j++){
                if (data[j].executed == 1 && fq_nmod_mpoly_is_zero(data[j].S_poly, ctx))
                    if (cmpPair(&sp, &(data[j].pPair)) == 0){
                        g_array_remove_index(P, i);
                        flag = 1;
                        break;
                    }

            }
            if (flag == 1) continue;
            i++;
        }
        return;
    }

//----------------------------------------------------
    _P = g_array_new(FALSE, FALSE, sizeof(Pair));
    rem_items = g_array_new(FALSE, FALSE, sizeof(ulong));
    fq_nmod_mpoly_init(lcm, ctx);
    fq_nmod_mpoly_init(div, ctx);
    fq_nmod_mpoly_init(lt_h, ctx);
    fq_nmod_mpoly_init(L, ctx);
    fq_nmod_mpoly_init(lt_f, ctx);
    fq_nmod_mpoly_init(gcd, ctx);

    LT(lt_h, h, ctx);
//----------------------------------------------------
    for(i = 0; i < t; i++){
        Pair sp = {i, t};
        g_array_append_val(_P, sp);
    }

    i = 0;
    
    while(i < P->len){
        Pair sp = g_array_index(P, Pair, i);
        flag = 0;
        for(j=0; j<data_len; j++){
            // print_poly("LOL", data[j].S_poly, NULL, ctx);
            if (data[j].executed != 1 || fq_nmod_mpoly_is_zero(data[j].S_poly, ctx) != 1) continue;
            if (cmpPair(&sp, &(data[j].pPair)) == 0){
                g_array_remove_index(P, i);
                flag = 1;
                break;
            }
        }

        

        if (flag == 1) continue;

        f = g_array_index(F, Polynom, sp.first);
        g = g_array_index(F, Polynom, sp.second);
        LCM(L, f, g, ctx);

        if (fq_nmod_mpoly_divides(div, L, lt_h, ctx) == 1){
            LCM(lcm, h, f, ctx);
            if (fq_nmod_mpoly_equal(lcm, L, ctx) == 0){
                LCM(lcm, h, g, ctx);
                if (fq_nmod_mpoly_equal(lcm, L, ctx) == 0){
                    g_array_remove_index(P, i);
                    i--;
                }
            }
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, Pair, i).first);
        LCM(lcm, f, h, ctx);
        j = 0;
        while(j < _P->len){
            if (i != j){
                
                g = g_array_index(F, Polynom, g_array_index(_P, Pair, j).first);
                LCM(L, g, h, ctx);
                if (fq_nmod_mpoly_divides(div, L, lcm, ctx) == 1){
                    g_array_remove_index(_P, j);

                    if (j < i)
                        i--;
                }
            }
            j++;
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        f = g_array_index(F, Polynom, g_array_index(_P, Pair, i).first);
        LT(lt_f, f, ctx);

        fq_nmod_mpoly_gcd(gcd, lt_f, lt_h, ctx);
        if (fq_nmod_mpoly_is_one(gcd, ctx) == 1){
            g_array_remove_index(_P, i);
        }
        else
            i++;
    }

    if (t == F->len){
        Polynom new_poly = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(new_poly, ctx);
        fq_nmod_mpoly_set(new_poly, h, ctx);

        g_array_append_val(F, new_poly);
    }

    pspair = (Pair*)_P->data;
    for(i = 0; i < _P->len; i++){
        Pair sp = {pspair->first, t};
        g_array_append_val(P, sp);
        pspair++;
    }

//----------------------------------------------------
    g_array_free(_P, TRUE);
    g_array_free(rem_items, TRUE);
    fq_nmod_mpoly_clear(lcm, ctx);
    fq_nmod_mpoly_clear(div, ctx);
    fq_nmod_mpoly_clear(lt_h, ctx);
    fq_nmod_mpoly_clear(L, ctx);
    fq_nmod_mpoly_clear(lt_f, ctx);
    fq_nmod_mpoly_clear(gcd, ctx);
}

Buchberger_result threaded_buchberger_v2(const Basis basis, ulong t, ulong threads_count, PolynomRing ctx){
    GArray *F, *P;
    Polynom* hp;
    Basis S_polys;
    thread_buchberger_data_v2_t* data;
    pthread_t* threads;
    int i, j, k;
    ulong executed_threads;
    ulong* selected_pairs;
    Pair sp;
    Polynom f, g;
//--------------------------------------------------------------
    F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    P = g_array_new(FALSE, FALSE, sizeof(Pair));
    S_polys = init_empty_basis(threads_count, ctx);
    data = flint_calloc(threads_count, sizeof(thread_buchberger_data_v2_t));
    threads = flint_calloc(threads_count, sizeof(pthread_t));
    selected_pairs = flint_calloc(threads_count, sizeof(ulong));

    for(i = 0; i < t; i++){
        Polynom f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_set(f, basis[i], ctx);
        g_array_append_val(F, f);
    }

    for(i = 0; i < threads_count; i++){
        data[i].ctx = ctx;
        data[i].executed = 0;
        data[i].executed_threads = &executed_threads;
        data[i].F = F;
        data[i].pause = 1;
        data[i].running = 1;
        data[i].S_poly = S_polys[i];

        pthread_create(&threads[i], NULL, buchberger_task_v2, &(data[i]));
    }
//--------------------------------------------------------------
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        GMI_v2(F, P, *hp, i, ctx);
        hp++;
    }

    // print_basis((Basis)F->data, F->len, NULL, ctx);

    while(P->len > 0){
        // print_basis((Basis)F->data, F->len, NULL, ctx);

        // printf("Pairs:\n");
        // for(i = 0; i < P->len; i++){
        //     Pair p = g_array_index(P, Pair, i);
        //     printf("(%ld, %ld)\n", p.first, p.second);
        // }
        // printf("\n");

        find_min_v3(F, P, selected_pairs, threads_count, ctx);

        // printf("Selected pairs:\n");
        // for(i = 0; i < threads_count; i++){
        //     if (selected_pairs[i] == -1) continue;
        //     Pair p = g_array_index(P, Pair, selected_pairs[i]);
        //     printf("(%ld, %ld)\n", p.first, p.second);
        // }
        // printf("\n");

        j = 0;
        for(i = 0; i < threads_count; i++){
            k = selected_pairs[i];
            data[i].executed = 0;
            if (k != -1){
                sp = g_array_index(P, Pair, k);
                f = g_array_index(F, Polynom, sp.first);
                g = g_array_index(F, Polynom, sp.second);
                data[i].pair = k;
                data[i].pPair = sp;
                S(S_polys[i], f, g, ctx);
                j++;
            }
        }
            
        if (j < threads_count){
            for(j; j < threads_count; j++)
                fq_nmod_mpoly_zero(S_polys[j], ctx);
        }

        // printf("S polys:\n");
        // for(i = 0; i < threads_count; i++){
        //     fq_nmod_mpoly_print_pretty(data[i].S_poly, NULL, ctx);
        //     printf("\n");
        // }
        // printf("\n");

        // print_basis((Basis)F->data, F->len, NULL, ctx);

        executed_threads = 0;
        for(i = 0; i < threads_count; i++){
            if (selected_pairs[i] == -1) continue;
            data[i].pause = 0;
        }

        // fq_nmod_mpoly_t reminder;
        // Basis Q;
        // fq_nmod_mpoly_init(reminder, ctx);
        // Q = init_empty_basis(F->len, ctx);
        // fq_nmod_mpoly_divrem_ideal(Q, reminder, S_polys[0], (Basis)F->data, F->len, ctx);
        // printf("remider[0]=");
        // fq_nmod_mpoly_print_pretty(reminder, NULL, ctx);
        // printf("\n");
        // fq_nmod_mpoly_clear(reminder, ctx);
        // free_basis(Q, F->len, ctx);

        int flag;
        while(executed_threads < NO_OF_IRRED){
            flag = 1;
            for(i = 0; i < threads_count; i++){
                if (data[i].pause == 0){
                    flag = 0;
                    break;
                }
            }
            if (flag == 1) break;
        }

        // printf("reminders:\n");
        // for(i = 0; i < threads_count; i++){
        //     fq_nmod_mpoly_print_pretty(data[i].S_poly, NULL, ctx);
        //     printf(" executed=%d\n", data[i].executed);
        // }
        // printf("\n");


        k = find_min_reduction_v2(data, threads_count, ctx);
        // printf("Selected reminder: ");
        // if (k != -1)
        //     fq_nmod_mpoly_print_pretty(data[k].S_poly, NULL, ctx);
        // else
        //     printf("k=-1");
        // printf("\n");
        

        // if (k != -1){
        //     for(i = 0; i < threads_count; i++){
        //         if (data[i].pair == k) continue;
        //         if (data[i].executed == 1 && data[i].pair > k) data[i].pair--;
        //     }
        //     g_array_remove_index(P, selected_pairs[k]);
        //     GMI_v3(F, P, S_polys[k], F->len, data, threads_count, ctx);
        // } else {

        // }

        if (k!=-1){
            g_array_remove_index(P, selected_pairs[k]);
            // printf("LOL\n");
            GMI_v3(F, P, S_polys[k], F->len, data, threads_count, ctx);
        } else {
            GMI_v3(F, P, NULL, F->len, data, threads_count, ctx);
        }



        // g_array_sort(P, cmpPair);
        // // printf("LOL\n");

        // for(i = 0; i < threads_count; i++){
        //     if(data[i].executed == 1 && fq_nmod_mpoly_is_zero(data[i].S_poly, ctx) == 1){
        //         if (g_array_binary_search(P, &(data[i].pPair), cmpPair, &k) == 1)
        //             g_array_remove_index(P, k);
        //     }
        // }

        // for(i = 0; i < threads_count; i++){
        //     if (data[i].executed == 1 && fq_nmod_mpoly_is_zero(data[i].S_poly, ctx) == 1){
        //         for(j = 0; j < threads_count; j++){
        //             if (i == j) continue;
        //             if (data[j].pair > data[i].pair) data[j].pair--;
        //         }
        //         Pair p = g_array_index(P, Pair, data[i].pair);
        //         // printf("remove pair (%ld, %ld) because S=", p.first, p.second);
        //         // fq_nmod_mpoly_print_pretty(data[i].S_poly, NULL, ctx);
        //         // printf("\n");
        //         printf("index: %ld %d; (%ld, %ld) (%ld, %ld)\n", data[i].pair, P->len, data[i].pPair.first, data[i].pPair.second, p.first, p.second);
        //         g_array_remove_index(P, data[i].pair);
        //     }
        // }
        
        // printf("\n");

        // sleep(5);

        // break;
    }
//--------------------------------------------------------------
    for(i = 0; i < threads_count; i++){
        data[i].running = 0;
        pthread_join(threads[i], NULL);
    }
    Basis res = from_garray(F);
    Buchberger_result resres = {res, F->len};
    g_array_free(F, TRUE);
    g_array_free(P, TRUE);
    free_basis(S_polys, threads_count, ctx);
    flint_free(data);
    flint_free(threads);
    flint_free(selected_pairs);

    return resres;
}
