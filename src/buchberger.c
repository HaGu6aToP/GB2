#include "headers.h"

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

// 1 - если критерий выполняется, 0 - не выполняется
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
            // log_S(S_polynom, f, g, ctx);
            fq_nmod_mpoly_divrem_ideal(Q, S_mod_G, S_polynom, basis, len, ctx);
            if (fq_nmod_mpoly_is_zero(S_mod_G, ctx) == 0){
                
                // printf("\n\ni----------------------------is_groebner_basis----------------------------\n");
                // printf("f: ");
                // fq_nmod_mpoly_print_pretty(f, NULL, ctx);
                // printf("\ng: ");
                // fq_nmod_mpoly_print_pretty(g, NULL, ctx);
                // printf("\nS(f, g): ");
                // fq_nmod_mpoly_print_pretty(S_polynom, NULL, ctx);
                // printf("\nS mod G: ");
                // fq_nmod_mpoly_print_pretty(S_mod_G, NULL, ctx);
                // printf("\n\n");

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

        flag1 = fq_nmod_mpoly_divides(div, L, lt_h, ctx);

        LCM(lcm, h, f, ctx);
        flag2 = fq_nmod_mpoly_equal(lcm, L, ctx);

        LCM(lcm, h, g, ctx);
        flag3 = fq_nmod_mpoly_equal(lcm, L, ctx);

        if((flag1 == 1) && (flag2 == 0) && (flag3 == 0)){
            free_SPair(&sp, ctx);
            g_array_remove_index(P, i);
            i--;
        }
        i++;
    }

    i = 0;
    while(i < _P->len){
        
        j = 0;
        while(j < _P->len){
            if (i != j){
                f = g_array_index(F, Polynom, g_array_index(_P, SPair, i).first);
                g = g_array_index(F, Polynom, g_array_index(_P, SPair, j).first);

                LCM(lcm, f, h, ctx);
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