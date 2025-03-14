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
        flag2 = crit(G, B, i, j, f, g, lcm, ctx);

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

