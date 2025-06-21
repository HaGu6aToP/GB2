#include "f4.h"
#include <unistd.h>
#include "flint/fq_nmod_mat.h"

#include "tools.h"
#include "basis_tools.h"
#include "sparse_matrix.h"

#define __DEBUG_F4 0

// #include <filesystem>
// #include <fstream>
// #include <iostream>

// #include "../SparseRREF/argparse.hpp"
// #include "../SparseRREF/sparse_mat.h"

#define HM(res, f, ctx) fq_nmod_mpoly_get_term_monomial(res, f, 0, ctx)
#define HT(res, f, ctx) fq_nmod_mpoly_get_term(res, f, 0, ctx)
#define HC(res, f, ctx) fq_nmod_mpoly_get_term_coeff_fq_nmod(res, f, 0, ctx)

#define init_poly(f, ctx) fq_nmod_mpoly_init(f, ctx)
#define clear_poly(f, ctx) fq_nmod_mpoly_clear(f, ctx)
#define set_poly(res, f, ctx) fq_nmod_mpoly_set(res, f, ctx)

void g_array_pop(GArray* p){
    g_array_remove_index(p, p->len-1);
}

void lcm(Polynom res, const Polynom f, const Polynom g, const PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    ulong exp_p1[nvars];
    ulong exp_p2[nvars];
    ulong exp_monom[nvars];


    fq_nmod_mpoly_get_term_exp_ui(exp_p1, f, 0, ctx);
    fq_nmod_mpoly_get_term_exp_ui(exp_p2, g, 0, ctx);

    for (int i = 0; i < nvars; i++)
        exp_monom[i] = max(exp_p1[i], exp_p2[i]);

    fq_nmod_mpoly_one(res, ctx);
    fq_nmod_mpoly_set_term_exp_ui(res, 0, exp_monom, ctx);
}

void spol(Polynom res, const Polynom f, const Polynom g, const Field field, const PolynomRing ring){
    fq_nmod_t hc_f, hc_g;
    fq_nmod_mpoly_t lcm_poly, hm_f, hm_g, temp, temp2;
    fq_nmod_init(hc_f, field);
    fq_nmod_init(hc_g, field);
    init_poly(lcm_poly, ring);
    init_poly(hm_f, ring);
    init_poly(hm_g, ring);
    init_poly(temp, ring);
    init_poly(temp2, ring);
    
    HC(hc_f, f, ring);
    HC(hc_g, g, ring);
    HM(hm_f, f, ring);
    HM(hm_g, g, ring);
    lcm(lcm_poly, f, g, ring);
    // print_poly("lcm", lcm_poly, NULL, ring);
    // print_poly("hm_f", hm_f, NULL, ring);
    // print_poly("hm_g", hm_g, NULL, ring);

    fq_nmod_mpoly_div(temp, lcm_poly, hm_f, ring); //lcm(f, g)/HM(f)
    fq_nmod_mpoly_scalar_mul_fq_nmod(res, temp, hc_g, ring); //HC(g)*lcm(f, g)/HM(f) 
    fq_nmod_mpoly_mul(temp, res, f, ring); //HC(g)*lcm(f, g)/HM(f) * f

    // fq_nmod_mpoly_print_pretty(temp, NULL, ring);
    // printf("\n");

    fq_nmod_mpoly_div(temp2, lcm_poly, hm_g, ring); 
    fq_nmod_mpoly_scalar_mul_fq_nmod(res, temp2, hc_f, ring); 
    fq_nmod_mpoly_mul(temp2, res, g, ring); 

    // fq_nmod_mpoly_print_pretty(temp2, NULL, ring);
    // printf("\n");

    fq_nmod_mpoly_sub(res, temp, temp2, ring); 

    fq_nmod_clear(hc_f, field);
    fq_nmod_clear(hc_g, field);
    clear_poly(lcm_poly, ring);
    clear_poly(hm_f, ring);
    clear_poly(hm_g, ring);
    clear_poly(temp, ring);
    clear_poly(temp2, ring);
}

void spol_old(Polynom res, const Polynom p1, const Polynom p2, const PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    fq_nmod_mpoly_t lcm_poly;
    fq_nmod_mpoly_init(lcm_poly, ctx);

    fq_nmod_mpoly_one(lcm_poly, ctx);
    lcm(lcm_poly, p1, p2, ctx);

    fq_nmod_mpoly_t leading_monom_p1, leading_monom_p2, A;
    fq_nmod_mpoly_init(leading_monom_p1, ctx);
    fq_nmod_mpoly_init(leading_monom_p2, ctx);
    fq_nmod_mpoly_init(A, ctx);

    fq_nmod_mpoly_get_term(leading_monom_p1, p1, 0, ctx);
    fq_nmod_mpoly_get_term(leading_monom_p2, p2, 0, ctx);

    fq_nmod_mpoly_div(A, lcm_poly, leading_monom_p1, ctx); // LCM(p1, p2) / LT(p1)
    fq_nmod_mpoly_mul(leading_monom_p1, A, p1, ctx); // LCM(p1, p2) / LT(f) * p1

    fq_nmod_mpoly_div(A, lcm_poly, leading_monom_p2, ctx); //LCM(p1, p2) / LT(p1)
    fq_nmod_mpoly_mul(leading_monom_p2, A, p2, ctx); //LCM(p1, p2) / LT(f) * p2

    fq_nmod_mpoly_sub(res, leading_monom_p1, leading_monom_p2, ctx); // S

    fq_nmod_mpoly_clear(lcm_poly, ctx);
    fq_nmod_mpoly_clear(leading_monom_p1, ctx);
    fq_nmod_mpoly_clear(leading_monom_p2, ctx);
    fq_nmod_mpoly_clear(A, ctx);
}

ulong deg(const Polynom f, PolynomRing ctx){
    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    ulong exp[nvars];
    fq_nmod_mpoly_get_term_exp_ui(exp, f, 0, ctx);
    return sum(exp, nvars);
}

void print_F4Pair(const F4Pair* p, const PolynomRing ctx){
    printf("{ ");
    printf("lcm=");
    fq_nmod_mpoly_print_pretty(p->lcm, NULL, ctx);
    printf(", t_f=");
    fq_nmod_mpoly_print_pretty(p->t_f, NULL, ctx);
    printf(", f=");
    fq_nmod_mpoly_print_pretty(p->f, NULL, ctx);
    printf(", t_g=");
    fq_nmod_mpoly_print_pretty(p->t_g, NULL, ctx);
    printf(", g=");
    fq_nmod_mpoly_print_pretty(p->g, NULL, ctx);
    printf(" }");
}

void print_F4PairProjection(const F4PairProjection* p, PolynomRing ctx){
    printf("{ t=");
    fq_nmod_mpoly_print_pretty(p->t, NULL, ctx);
    printf(", f=");
    fq_nmod_mpoly_print_pretty(p->f, NULL, ctx);
    printf(" }\n");
}

void init_F4Pair(F4Pair* p, const Polynom f, const Polynom g, const PolynomRing ctx){
    p->f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    p->g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    p->lcm = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    p->t_f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    p->t_g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));

    init_poly(p->f, ctx);
    init_poly(p->g, ctx);
    init_poly(p->lcm, ctx);
    init_poly(p->t_f, ctx);
    init_poly(p->t_g, ctx);

    set_poly(p->f, f, ctx);
    set_poly(p->g, g, ctx);
    lcm(p->lcm, f, g, ctx);

    ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    ulong exp_hm_f[nvars];
    ulong exp_hm_g[nvars];
    ulong exp_hm_lcm[nvars];

    fq_nmod_mpoly_get_term_exp_ui(exp_hm_f, f, 0, ctx);
    fq_nmod_mpoly_get_term_exp_ui(exp_hm_g, g, 0, ctx);
    fq_nmod_mpoly_get_term_exp_ui(exp_hm_lcm, p->lcm, 0, ctx);

    for(int i = 0; i < nvars; i++){
        exp_hm_f[i] = exp_hm_lcm[i] - exp_hm_f[i];
        exp_hm_g[i] = exp_hm_lcm[i] - exp_hm_g[i];
    }

    fq_nmod_mpoly_one(p->t_f, ctx);
    fq_nmod_mpoly_one(p->t_g, ctx);

    fq_nmod_mpoly_set_term_exp_ui(p->t_f, 0, exp_hm_f, ctx);
    fq_nmod_mpoly_set_term_exp_ui(p->t_g, 0, exp_hm_g, ctx);
}

void free_F4Pair(F4Pair* p, const PolynomRing ctx){
    clear_poly(p->f, ctx);
    clear_poly(p->g, ctx);
    clear_poly(p->lcm, ctx);
    clear_poly(p->t_f, ctx);
    clear_poly(p->t_g, ctx);

    flint_free(p->f);
    flint_free(p->g);
    flint_free(p->lcm);
    flint_free(p->t_f);
    flint_free(p->t_g);
}

void free_F4Pair_lst(GArray* P, const PolynomRing ctx){
    F4Pair f4p;

    while(P->len > 0){
        f4p = g_array_index(P, F4Pair, P->len-1);
        free_F4Pair(&f4p, ctx);
        g_array_remove_index(P, P->len-1);
    }

    g_array_free(P, TRUE);
}

void print_F4Pair_lst(const GArray* P, const PolynomRing ctx){
    if (P->len == 0) return;
    for(ulong i = 0; i < P->len-1; i++){
        print_F4Pair(&g_array_index(P, F4Pair, i), ctx);
        printf("\n");
    }
    print_F4Pair(&g_array_index(P, F4Pair, P->len-1), ctx);
}

ulong find_min_deg_in_F4Pairs(const GArray* P, const PolynomRing ctx){
    if (P->len == 0)
        return 0;

    

    F4Pair* pf4p = (F4Pair*)P->data;
    ulong d = deg(pf4p->lcm, ctx);
    ulong a = 0;
    pf4p++;

    for(ulong i = 1; i < P->len; i++){
        a = deg(pf4p->lcm, ctx);
        if (a < d) d = a;
        pf4p++;
    }

    return d;
}

void F4_select(GArray* Pd, GArray* P, const PolynomRing ctx){
    // if (Pd != NULL) g_array_free(Pd, TRUE);
    // Pd = g_array_new(FALSE, FALSE, sizeof(F4Pair));
    
    if (P == NULL || P->len == 0)
        return;

    ulong d = find_min_deg_in_F4Pairs(P, ctx);
    F4Pair f4p;
    ulong i = 0;

    while(i < P->len){
        f4p = g_array_index(P, F4Pair, i);
        if (deg(f4p.lcm, ctx) == d){
            g_array_remove_index(P, i);
            g_array_append_val(Pd, f4p);
            continue;
        }
        i++;
    }
}

// preprocessing с массивами
void old_old_preprocessing(GArray* F, GArray* Pd, const GArray* G, const PolynomRing ctx){
    F4Pair f4p, *pf4p;
    Polynom new_poly, f;
    Polynom* hp;
    fq_nmod_mpoly_t m, div;
    GArray* done;
    GArray* sub;
    ulong i, j, k;
//-------------------------------------------------------
    done = __calloc_poly_lst();
    sub = __calloc_poly_lst();
    init_poly(m, ctx);
    init_poly(div, ctx);

    pf4p = (F4Pair*)Pd->data;
    for(i = 0; i < Pd->len; i++){
        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, pf4p[i].t_f, pf4p[i].f, ctx);
        g_array_append_val(F, new_poly);

        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, pf4p[i].t_g, pf4p[i].g, ctx);
        g_array_append_val(F, new_poly);

        free_F4Pair(&pf4p[i], ctx);
    }
    g_array_remove_range(Pd, 0, Pd->len);

    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        if (i % 2 == 0){
            HM(m, *hp, ctx);
            if (is_poly_in_lst(done, m, ctx) == 0){
                f = __calloc_poly();
                init_poly(f, ctx);
                set_poly(f, m, ctx);
                g_array_append_val(done, f);
            }
        }

        for(j = 1; j < fq_nmod_mpoly_length(*hp, ctx); j++){
            fq_nmod_mpoly_get_term_monomial(m, *hp, j, ctx);
            if (is_poly_in_lst(sub, m, ctx) == 0){
                f = __calloc_poly();
                init_poly(f, ctx);
                set_poly(f, m, ctx);
                g_array_append_val(sub, f);
            }
        }

        hp++;
    }
//-------------------------------------------------------
    #if __DEBIG
        printf("---------------------------------------preprocessing---------------------------------------\n");
    #endif
    
    // Добавляем новые полиномы до тех пор, пока для каждого монома из T(F) не найдется полином f из F для котрого он ведущий
    while(sub->len != 0){
        #if __DEBUG_F4
            printf("F:\n");
            print_poly_lst(F, ctx);
            printf("\n");

            printf("Done:\n");
            print_poly_lst(done, ctx);
            printf("\n");

            printf("HM(F)\\Done:\n");
            print_poly_lst(sub, ctx);
            printf("\n");
        #endif

        k = max_poly_in_lst(sub, ctx);
        #if __DEBUG_F4
            printf("k=%ld\n", k);
        #endif

        f = g_array_index(sub, Polynom, k);

        #if __DEBUG_F4
            printf("selected monom: ");
            fq_nmod_mpoly_print_pretty(f, NULL, ctx);
            printf("\n");
        #endif

        g_array_append_val(done, f);
        g_array_remove_index(sub, k);
        
        hp = (Polynom*)G->data;
        for(i = 0; i < G->len; i++){
            HM(m, *hp, ctx);
            if (fq_nmod_mpoly_divides(div, f, m, ctx) == 1){
                #if __DEBUG_F4
                    printf("selected monom - ");
                    fq_nmod_mpoly_print_pretty(f, NULL, ctx);
                    printf(" divides by HT(");
                    fq_nmod_mpoly_print_pretty(*hp, NULL, ctx);
                    printf(")=");
                    fq_nmod_mpoly_print_pretty(m, NULL, ctx);
                    printf("\n");
                    printf("div=");
                    fq_nmod_mpoly_print_pretty(div, NULL, ctx);
                    printf("\n");
                #endif

                new_poly = __calloc_poly();
                init_poly(new_poly, ctx);
                fq_nmod_mpoly_mul(new_poly, div, *hp, ctx);
                g_array_append_val(F, new_poly);

                #if __DEBUG_F4
                    printf("new poly - ");
                    fq_nmod_mpoly_print_pretty(new_poly, NULL, ctx);
                    printf("\n");
                #endif

                for(i = 1; i < fq_nmod_mpoly_length(new_poly, ctx); i++){
                    fq_nmod_mpoly_get_term_monomial(m, new_poly, i, ctx);
                    if (is_poly_in_lst(sub, m, ctx) == 0){
                        f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
                        init_poly(f, ctx);
                        fq_nmod_mpoly_get_term_monomial(f, new_poly, i, ctx);
                        g_array_append_val(sub, f);
                    }
                }
                break;
            }
            hp++;
        }

        // break;
    }
    #if __DEBUG_F4
        printf("---------------------------------------preprocessing---------------------------------------\n");
    #endif
//-------------------------------------------------------
    // g_array_free(sub, TRUE);
    free_poly_lst(sub, ctx);
    free_poly_lst(done, ctx);
    clear_poly(m, ctx);
    clear_poly(div, ctx);
}

// preprocessing с хеш-таблицами вместо массивов
void old_preprocessing(GArray* F, GArray* Pd, const GArray* G, const PolynomRing ctx){
    F4Pair f4p, *pf4p;
    Polynom new_poly, f;
    Polynom* hp;
    Polynom h;
    fq_nmod_mpoly_t m, div;
    GHashTable* done;
    GHashTable* sub;
    gpointer* gp;
    ulong i, j, k;
    char* str;
// -------------------------------------------------
    done = g_hash_table_new_full(g_str_hash, g_str_equal, simple_key_destroyer, NULL);
    sub =  g_hash_table_new_full(g_str_hash, g_str_equal, simple_key_destroyer, NULL);
    init_poly(m, ctx);
    init_poly(div, ctx);
    

    // Формирование множеств Left и Right
    pf4p = (F4Pair*)Pd->data;
    for(i = 0; i < Pd->len; i++){
        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, pf4p[i].t_f, pf4p[i].f, ctx);
        g_array_append_val(F, new_poly);

        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, pf4p[i].t_g, pf4p[i].g, ctx);
        g_array_append_val(F, new_poly);

        free_F4Pair(&pf4p[i], ctx);
    }
    g_array_remove_range(Pd, 0, Pd->len);

    // Формирование множества done и sub
    for(i = 0; i < F->len; i++){
        h = g_array_index(F, Polynom, i);
        HM(m, h, ctx);

        // Получаем строковое представление монома для вычисления хеша
        str = fq_nmod_mpoly_get_str_pretty(m, NULL, ctx);
        // printf("str = %s\n", str);
        // Проверяем есть ли уже этот моном, если нет добовляем
        if (g_hash_table_lookup(done, str) == NULL){
            // printf("inserting ");
            f = __calloc_poly();
            init_poly(f, ctx);
            set_poly(f, m, ctx);
            // printf("str = %s\n", str);
            g_hash_table_insert(done, str, f);

            // printf("done(%d): \n", g_hash_table_size(done));
            // print_hash_table(done, ctx);
            // printf("\n");
        }
        // free(str);

        for(j = 1; j < fq_nmod_mpoly_length(h, ctx); j++){
            fq_nmod_mpoly_get_term_monomial(m, h, j, ctx);

            str = fq_nmod_mpoly_get_str_pretty(m, NULL, ctx);

            if (g_hash_table_lookup(sub, str) == NULL){
                f = __calloc_poly();
                init_poly(f, ctx);
                set_poly(f, m, ctx);
                g_hash_table_insert(sub, str, f);
            }
            // free(str);
        }
    }

    #if __DEBUG_F4
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\n");
        printf("done: \n");
        print_hash_table(done, ctx);
        printf("\nsub: \n");
        print_hash_table(sub, ctx);
        printf("\n");
    #endif

    // Основной цикл
    while(g_hash_table_size(sub) > 0){


        // Находим наибольший моном и перемещаем его из sub в done
        f = max_poly_in_GHashtable(sub, ctx);

        #if __DEBUG_F4
            printf("Наибольший моном: ");
            fq_nmod_mpoly_print_pretty(f, NULL, ctx);
            printf("\n");
        #endif

        str = fq_nmod_mpoly_get_str_pretty(f, NULL, ctx);
        // if (g_hash_table_size(sub) == 1)
        //     g_hash_table_remove_all(sub);
        // else
            g_hash_table_remove(sub, str);
        g_hash_table_insert(done, str, f);

        #if __DEBUG_F4
            printf("done(%d): \n", g_hash_table_size(done));
            print_hash_table(done, ctx);
            printf("\nsub(%d): \n", g_hash_table_size(sub));
            print_hash_table(sub, ctx);
            printf("\n");
            sleep(3);
        #endif

        // Проверяем моном f на делимость некоторым LM(G[i])
        for (i = 0; i < G->len; i++){
            h = g_array_index(G, Polynom, i);
            HM(m, h, ctx);


            if (fq_nmod_mpoly_divides(div, f, m, ctx) == 1){
                #if __DEBUG_F4
                    printf("selected monom - ");
                    fq_nmod_mpoly_print_pretty(f, NULL, ctx);
                    printf(" divides by HT(");
                    fq_nmod_mpoly_print_pretty(h, NULL, ctx);
                    printf(")=");
                    fq_nmod_mpoly_print_pretty(m, NULL, ctx);
                    printf("\n");
                    printf("div=");
                    fq_nmod_mpoly_print_pretty(div, NULL, ctx);
                    printf("\n");
                #endif

                // Добавляем новый полином
                new_poly = __calloc_poly();
                init_poly(new_poly, ctx);
                fq_nmod_mpoly_mul(new_poly, div, h, ctx);
                g_array_append_val(F, new_poly);

                #if __DEBUG_F4
                    printf("new poly - ");
                    fq_nmod_mpoly_print_pretty(new_poly, NULL, ctx);
                    printf("\n");
                #endif

                // Обновляем множество sub
                for(j = 1; j < fq_nmod_mpoly_length(new_poly, ctx); j++){
                    fq_nmod_mpoly_get_term_monomial(m, new_poly, j, ctx);

                    str = fq_nmod_mpoly_get_str_pretty(m, NULL, ctx);

                    if (g_hash_table_lookup(sub, str) == NULL){
                        f = __calloc_poly();
                        init_poly(f, ctx);
                        set_poly(f, m, ctx);
                        g_hash_table_insert(sub, str, f);
                    }
                    // free(str);
                }
                break;

            }
        }

    }

    // Особождение данных
    clear_poly(m, ctx);
    clear_poly(div, ctx);


    GPtrArray* vals = g_hash_table_get_values_as_ptr_array(done);
    for(i = 0; i < vals->len; i++){
        clear_poly((Polynom)vals->pdata[i], ctx);
    }

    g_ptr_array_free(vals, TRUE);

    g_hash_table_destroy(done);
    g_hash_table_destroy(sub);

}

// preprocessing с новым ключом
void preprocessing(GArray* F, GArray* Pd, const GArray* G, const PolynomRing ctx){
    F4Pair f4p, *pf4p;
    Polynom new_poly, f;
    Polynom* hp;
    Polynom h;
    fq_nmod_mpoly_t m, div;
    GHashTable* done;
    GHashTable* sub;
    gpointer* gp;
    ulong i, j, k;
    ulong *key, *new_key, ulong_key;
// -------------------------------------------------
    done = g_hash_table_new_full(g_int64_hash, g_int64_equal, simple_key_destroyer, NULL);
    sub =  g_hash_table_new_full(g_int64_hash, g_int64_equal, simple_key_destroyer, NULL);
    init_poly(m, ctx);
    init_poly(div, ctx);
    

    // Формирование множеств Left и Right
    pf4p = (F4Pair*)Pd->data;
    for(i = 0; i < Pd->len; i++){
        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, pf4p[i].t_f, pf4p[i].f, ctx);
        g_array_append_val(F, new_poly);

        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, pf4p[i].t_g, pf4p[i].g, ctx);
        g_array_append_val(F, new_poly);

        free_F4Pair(&pf4p[i], ctx);
    }
    g_array_remove_range(Pd, 0, Pd->len);

    // Формирование множества done и sub
    for(i = 0; i < F->len; i++){
        h = g_array_index(F, Polynom, i);
        HM(m, h, ctx);

        // Получаем ключ
        ulong_key = monom_hash(m, ctx);
        
        // Проверяем есть ли уже этот моном, если нет добовляем
        if (g_hash_table_lookup(done, &ulong_key) == NULL){
            // printf("inserting ");
            f = __calloc_poly();
            init_poly(f, ctx);
            set_poly(f, m, ctx);

            key = malloc(sizeof(ulong));
            *key = ulong_key;

            // printf("str = %s\n", str);
            g_hash_table_insert(done, key, f);
            // printf("done(%d): \n", g_hash_table_size(done));
            // print_hash_table(done, ctx);
            // printf("\n");
        }

        for(j = 1; j < fq_nmod_mpoly_length(h, ctx); j++){
            fq_nmod_mpoly_get_term_monomial(m, h, j, ctx);

            ulong_key = monom_hash(m, ctx);

            if (g_hash_table_lookup(sub, &ulong_key) == NULL){
                f = __calloc_poly();
                init_poly(f, ctx);
                set_poly(f, m, ctx);

                key = malloc(sizeof(ulong));
                *key = ulong_key;

                g_hash_table_insert(sub, key, f);
            }
        }
    }

    #if __DEBUG_F4
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\n");
        printf("done: \n");
        print_hash_table(done, ctx);
        printf("\nsub: \n");
        print_hash_table(sub, ctx);
        printf("\n");
    #endif


    // Основной цикл
    while(g_hash_table_size(sub) > 0){


        // Находим наибольший моном и перемещаем его из sub в done
        f = max_poly_in_GHashtable(sub, ctx);

        #if __DEBUG_F4
            printf("Наибольший моном: ");
            fq_nmod_mpoly_print_pretty(f, NULL, ctx);
            printf("\n");
        #endif

        key = malloc(sizeof(ulong));
        *key = monom_hash(f, ctx);
        // if (g_hash_table_size(sub) == 1)
        //     g_hash_table_remove_all(sub);
        // else
            g_hash_table_remove(sub, key);
        g_hash_table_insert(done, key, f);

        #if __DEBUG_F4
            printf("done(%d): \n", g_hash_table_size(done));
            print_hash_table(done, ctx);
            printf("\nsub(%d): \n", g_hash_table_size(sub));
            print_hash_table(sub, ctx);
            printf("\n");
            sleep(3);
        #endif

        // Проверяем моном f на делимость некоторым LM(G[i])
        for (i = 0; i < G->len; i++){
            h = g_array_index(G, Polynom, i);
            HM(m, h, ctx);


            if (fq_nmod_mpoly_divides(div, f, m, ctx) == 1){
                #if __DEBUG_F4
                    printf("selected monom - ");
                    fq_nmod_mpoly_print_pretty(f, NULL, ctx);
                    printf(" divides by HT(");
                    fq_nmod_mpoly_print_pretty(h, NULL, ctx);
                    printf(")=");
                    fq_nmod_mpoly_print_pretty(m, NULL, ctx);
                    printf("\n");
                    printf("div=");
                    fq_nmod_mpoly_print_pretty(div, NULL, ctx);
                    printf("\n");
                #endif

                // Добавляем новый полином
                new_poly = __calloc_poly();
                init_poly(new_poly, ctx);
                fq_nmod_mpoly_mul(new_poly, div, h, ctx);
                g_array_append_val(F, new_poly);

                #if __DEBUG_F4
                    printf("new poly - ");
                    fq_nmod_mpoly_print_pretty(new_poly, NULL, ctx);
                    printf("\n");
                #endif

                // Обновляем множество sub
                for(j = 1; j < fq_nmod_mpoly_length(new_poly, ctx); j++){
                    fq_nmod_mpoly_get_term_monomial(m, new_poly, j, ctx);

                    
                    ulong_key = monom_hash(m, ctx);

                    if (g_hash_table_lookup(sub, &ulong_key) == NULL){
                        f = __calloc_poly();
                        init_poly(f, ctx);
                        set_poly(f, m, ctx);

                        key = malloc(sizeof(ulong));
                        *key = ulong_key;

                        g_hash_table_insert(sub, key, f);
                    }
                }
                break;

            }
        }

    }


    // Особождение данных
    clear_poly(m, ctx);
    clear_poly(div, ctx);


    GPtrArray* vals = g_hash_table_get_values_as_ptr_array(done);
    for(i = 0; i < vals->len; i++){
        clear_poly((Polynom)vals->pdata[i], ctx);
    }

    g_ptr_array_free(vals, TRUE);

    g_hash_table_destroy(done);
    g_hash_table_destroy(sub);

}


// LinBox reduce
void ref4(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){
    GArray* F_monoms;
    Polynom f;
    Polynom *hp, *hg;
    fq_nmod_mpoly_t m, sum, g;
    fq_nmod_mat_t M;
    fq_nmod_t x, y;
    ulong i, j, k;
    slong r, p_len, t;
//-------------------------------------------------------
    F_monoms = __calloc_poly_lst();
    monom_lst_from_poly_lst(F_monoms, F, ctx);
    init_poly(m, ctx);
    init_poly(sum, ctx);
    init_poly(g, ctx);
    p_len = MAX(F_monoms->len, F->len);

    poly_quick_sort(F_monoms, 0, F_monoms->len-1, 1, ctx);
//-------------------------------------------------------

    #if __DEBUG_F4_POLY_REDUCE
            printf("F:\n");
            print_poly_lst(F, ctx);
            printf("\nF_monoms:\n");
            print_poly_lst(F_monoms, ctx);
            printf("\n");
    #endif

    F4_linbox_poly_reduce(F_ref, F, F_monoms, field, ctx);

    #if __DEBUG_F4_POLY_REDUCE
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\nF_monoms:\n");
        print_poly_lst(F_monoms, ctx);
        printf("\n");

        printf("F_ref:\n");
        print_poly_lst(F_ref, ctx);
        printf("\n");
        // sleep(1000);
    #endif
    
//-------------------------------------------------------
    free_poly_lst(F_monoms, ctx);
    clear_poly(m, ctx);
    clear_poly(sum, ctx);
    clear_poly(g, ctx);
}

// GLIB reduce
void ref3(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){
    GArray* F_monoms;
    Polynom f;
    Polynom *hp, *hg;
    fq_nmod_mpoly_t m, sum, g;
    fq_nmod_mat_t M;
    fq_nmod_t x, y;
    ulong i, j, k;
    slong r, p_len, t;
//-------------------------------------------------------
    F_monoms = __calloc_poly_lst();
    monom_lst_from_poly_lst(F_monoms, F, ctx);
    init_poly(m, ctx);
    init_poly(sum, ctx);
    init_poly(g, ctx);
    p_len = MAX(F_monoms->len, F->len);

    poly_quick_sort(F_monoms, 0, F_monoms->len-1, 1, ctx);
//-------------------------------------------------------

    #if __DEBUG_F4_POLY_REDUCE
            printf("F:\n");
            print_poly_lst(F, ctx);
            printf("\nF_monoms:\n");
            print_poly_lst(F_monoms, ctx);
            printf("\n");
    #endif

    F4_poly_reduce(F_ref, F, F_monoms, field, ctx);

    #if __DEBUG_F4_POLY_REDUCE
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\nF_monoms:\n");
        print_poly_lst(F_monoms, ctx);
        printf("\n");

        printf("F_ref:\n");
        print_poly_lst(F_ref, ctx);
        printf("\n");
        // sleep(1000);
    #endif
    
//-------------------------------------------------------
    free_poly_lst(F_monoms, ctx);
    clear_poly(m, ctx);
    clear_poly(sum, ctx);
    clear_poly(g, ctx);
}

void ref2(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){
    GArray* F_monoms;
    Polynom f;
    Polynom *hp, *hg;
    fq_nmod_mpoly_t m, sum, g;
    fq_nmod_mat_t M;
    fq_nmod_t x, y;
    ulong i, j, k;
    slong* p;
    slong r, p_len, t;
//-------------------------------------------------------
    F_monoms = __calloc_poly_lst();
    monom_lst_from_poly_lst(F_monoms, F, ctx);
    init_poly(m, ctx);
    init_poly(sum, ctx);
    init_poly(g, ctx);
    p_len = MAX(F_monoms->len, F->len);
    p = flint_calloc(p_len, sizeof(slong));
    // for(i = 0; i < p_len; i++)
    //     p[i] = i;

    poly_quick_sort(F_monoms, 0, F_monoms->len-1, 1, ctx);

    // printf("%d %d", F->len, F_monoms->len);

    fq_nmod_mat_init(M, F->len, F_monoms->len, field);
    fq_nmod_init(x, field);
    fq_nmod_init(y, field);

        // Формируем матрицу 
    // hp = (Polynom*)F->data;
    // for(i = 0; i < F->len; i++){
    //     for(j = 0; j < fq_nmod_mpoly_length(*hp, ctx); j++){
    //         fq_nmod_mpoly_get_term_monomial(m, *hp, j, ctx);
    //         fq_nmod_mpoly_get_term_coeff_fq_nmod(x, *hp, j, ctx);
            
    //         hg = (Polynom*)F_monoms->data;
    //         for(k = 0; k < F_monoms->len; k++){
    //             if (fq_nmod_mpoly_equal(m, *hg, ctx) == 1){
    //                 fq_nmod_mat_entry_set(M, i, k, x, field);
    //                 break;
    //             }
    //             hg++;
    //         }
    //     }
    //     hp++;
    // }

    // printf("M:\n");
    // fq_nmod_mat_print_pretty(M, field);
    // printf("\n");

    // r = fq_nmod_mat_lu_classical(p, M, 0, field);

    // printf("%d %d\n", F->len, F_monoms->len);
    // printf("M LU rank=%ld:\n", r);
    // fq_nmod_mat_print_pretty(M, field);
    // printf("\n");
    // // flint_free(p);

    // // printf("LOL %ld\n", p_len);
    // for(i = 0; i < p_len; i++)
    //     printf("%ld ", p[i]); 

    // slong tt;

    // // sleep(10);

    //    // Получаем редуцированные полиномы
    // for(i = 0; i < r; i++){
    //     f = __calloc_poly();
    //     init_poly(f, ctx);

    //     for(j=i; j < F_monoms->len; j++){
    //         if (fq_nmod_is_zero(fq_nmod_mat_entry(M, i, j), field) == 1) continue;

    //         if (F_monoms->len >= F->len && !(p[p_len - 1] == 0 && p[p_len - 2] == 0)) tt = p[j]; 
    //         else tt = j;

    //         // printf("k=%ld\n", tt);

    //         fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, tt), fq_nmod_mat_entry(M, i, j), ctx);
    //         set_poly(sum, f, ctx);
    //         fq_nmod_mpoly_add(f, sum, m, ctx);
    //     }

    //     // g_array_append_val(F_ref, f);
    //     fq_nmod_mpoly_print_pretty(f, NULL, ctx);
    //     printf("\n");
    //     clear_poly(f, ctx);
    //     flint_free(f);
    // }

//-------------------------------------------------------

    #if __DEBUG_F4
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\n");
    #endif
    

    sparse_matrix_t sparse_M;
    sparse_matrix_init(sparse_M, F->len, F_monoms->len, field);


    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        for(j = 0; j < fq_nmod_mpoly_length(*hp, ctx); j++){
            fq_nmod_mpoly_get_term_monomial(m, *hp, j, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(x, *hp, j, ctx);
            
            hg = (Polynom*)F_monoms->data;
            for(k = 0; k < F_monoms->len; k++){
                if (fq_nmod_mpoly_equal(m, *hg, ctx) == 1){
                    // fq_nmod_mat_entry_set(M, i, k, x, field);
                    sparse_matrix_add_elem_fq_nmod(sparse_M, i, k, x);
                    break;
                }
                hg++;
            }
        }
        hp++;
    }

    #if __DEBUG_F4
        sparse_matrix_print_pretty(sparse_M);
        printf("\n");
    #endif

    ulong not_null_elems = 0;
    for(i = 0; i < F->len; i++)
        not_null_elems += fq_nmod_mpoly_length(g_array_index(F, Polynom, i), ctx);

    #if __DEBUG_F4
        printf("LU (%d x %d) : %ld\n", F->len, F_monoms->len, not_null_elems);
    #endif

    r = sparse_matrix_gauss_rref(sparse_M);

    #if __DEBUG_F4
        printf("LU completed\n");

        printf("ref: r=%ld\n", r);
        sparse_matrix_print_pretty(sparse_M);
        printf("\n");

        sparse_matrix_print_info(sparse_M);
        printf("\n");

        printf("F_monoms:\n");
        print_poly_lst(F_monoms, ctx);
        printf("\n");
    #endif

    sparse_matrix_pair smp;
    

    for(i = 0; i < r; i++){
        f = __calloc_poly();
        init_poly(f, ctx);
        
        #if __DEBUG_F4
            printf("els = %ld\n", sparse_matrix_els_in_line(sparse_M, i));
        #endif

        for(j = 0; j < sparse_matrix_els_in_line(sparse_M, i); j++){
            // sparse_matrix_pair* smp = sparse_matrix_true_entry_fq_nmod(sparse_M, i, j);
            sparse_matrix_true_entry_fq_nmod(&smp, sparse_M, i, j);
            fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, smp.k), smp.val, ctx);
            // if (fq_nmod_is_zero(fq_nmod_mat_entry(M, i, j), field) == 1) continue;

            // if (F_monoms->len >= F->len && p[p_len - 1] != -1) tt = p[j]; 
            // else tt = j;

            // // printf("k=%ld\n", tt);

            // fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, tt), fq_nmod_mat_entry(M, i, j), ctx);
            set_poly(sum, f, ctx);
            fq_nmod_mpoly_add(f, sum, m, ctx);
        }

        g_array_append_val(F_ref, f);
        #if __DEBUG_F4
            fq_nmod_mpoly_print_pretty(f, NULL, ctx);
            printf("\n");
        #endif
        // clear_poly(f, ctx);
        // flint_free(f);
    }
    
    
//-------------------------------------------------------
    free_poly_lst(F_monoms, ctx);
    clear_poly(m, ctx);
    clear_poly(sum, ctx);
    clear_poly(g, ctx);
    fq_nmod_clear(x, field);
    fq_nmod_clear(y, field);
    fq_nmod_mat_clear(M, field);
    flint_free(p);
    sparse_matrix_clear(sparse_M);
}

// FLINT LU reduce
void ref(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){
    GArray* F_monoms;
    Polynom f;
    Polynom *hp, *hg;
    fq_nmod_mpoly_t m, sum, g;
    fq_nmod_mat_t M;
    fq_nmod_t x, y;
    ulong i, j, k;
    slong* p;
    slong r, p_len, t;
//-------------------------------------------------------
    F_monoms = __calloc_poly_lst();
    monom_lst_from_poly_lst(F_monoms, F, ctx);
    init_poly(m, ctx);
    init_poly(sum, ctx);
    init_poly(g, ctx);
    p_len = MAX(F_monoms->len, F->len);
    p = flint_calloc(p_len, sizeof(slong));
    // for(i = 0; i < p_len; i++)
    //     p[i] = i;

    poly_quick_sort(F_monoms, 0, F_monoms->len-1, 1, ctx);

    fq_nmod_mat_init(M, F->len, F_monoms->len, field);
    fq_nmod_init(x, field);
    fq_nmod_init(y, field);


    // Формируем матрицу 
    hp = (Polynom*)F->data;
    for(i = 0; i < F->len; i++){
        for(j = 0; j < fq_nmod_mpoly_length(*hp, ctx); j++){
            fq_nmod_mpoly_get_term_monomial(m, *hp, j, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(x, *hp, j, ctx);
            
            hg = (Polynom*)F_monoms->data;
            for(k = 0; k < F_monoms->len; k++){
                if (fq_nmod_mpoly_equal(m, *hg, ctx) == 1){
                    fq_nmod_mat_entry_set(M, i, k, x, field);
                    break;
                }
                hg++;
            }
        }
        hp++;
    }
    // printf("%d, %d\n", F->len, F_monoms->len);
//-------------------------------------------------------
    #if __DEBUG_F4
        printf("---------------------------------------ref---------------------------------------\n");
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\n");
        printf("F_monoms:\n");
        print_poly_lst(F_monoms, ctx);
        printf("\n");
        printf("%d %d\n", F_monoms->len, fq_nmod_mat_ncols(M, field));
        printf("Columns=%d, Lines=%d\n", F_monoms->len, F->len);
        printf("M:\n");
        fq_nmod_mat_print_pretty(M, field);
        printf("\n");
    
        printf("columns ordering: \n");
        for(i = 0; i < p_len; i++)
            printf("%ld ", p[i]);
        printf("\n");
    #endif

    // Приводим к верхне треугольней форме с помощью LU разложения
    ulong not_null_elems = 0;
    for(i = 0; i < F->len; i++)
        not_null_elems += fq_nmod_mpoly_length(g_array_index(F, Polynom, i), ctx);

    #if __DEBUG_F4
        printf("LU (%d x %d) : %ld\n", F->len, F_monoms->len, not_null_elems);
    #endif

    p[p_len - 1] = -1;
    // FILE* file = fopen("res.txt", "w");
    // fq_nmod_mat_fprint_pretty(file, M, field);
    // fclose(file);
    r = fq_nmod_mat_lu_classical(p, M, 0, field);

    #if __DEBUG_F4
        printf("LU complete\n");
        printf("new columns ordering: \n");
        for(i = 0; i < p_len; i++)
            printf("%ld ", p[i]);
        printf("\n");

        printf("\n");
        printf("%d %d\n", F_monoms->len, fq_nmod_mat_ncols(M, field));
        printf("M LU (%d, %d), %ld:\n", F->len, F_monoms->len, r);
        fq_nmod_mat_print_pretty(M, field);
        printf("\n");
        printf("rank=%ld\n", r);

        printf("\n");
        printf("M LU:\n");
        fq_nmod_mat_print_pretty(M, field);
        printf("\n");
    #endif

    // Получаем редуцированные полиномы
    for(i = 0; i < r; i++){
        f = __calloc_poly();
        init_poly(f, ctx);

        for(j=i; j < F_monoms->len; j++){
            if (fq_nmod_is_zero(fq_nmod_mat_entry(M, i, j), field) == 1) continue;

            if (F_monoms->len >= F->len && p[p_len - 1] != -1) t = p[j]; 
            else t = j;

            fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, t), fq_nmod_mat_entry(M, i, j), ctx);
            set_poly(sum, f, ctx);
            fq_nmod_mpoly_add(f, sum, m, ctx);
        }

        g_array_append_val(F_ref, f);
        #if __DEBUG_F4
            fq_nmod_mpoly_print_pretty(f, NULL, ctx);
            printf("\n");
        #endif
    }

    #if __DEBUG_F4
        printf("---------------------------------------ref-end---------------------------------------\n");
    #endif
//-------------------------------------------------------
    free_poly_lst(F_monoms, ctx);
    clear_poly(m, ctx);
    clear_poly(sum, ctx);
    clear_poly(g, ctx);
    fq_nmod_clear(x, field);
    fq_nmod_clear(y, field);
    fq_nmod_mat_clear(M, field);
    flint_free(p);
}

void reduction(GArray* F_, GArray* Pd, const GArray* G, const Field field, const PolynomRing ctx){
    GArray* F;
    GArray* F_ref;
    Polynom *hg;
    Polynom h;
    fq_nmod_mpoly_t f, g;
    ulong i, j;
    int flag;
//-------------------------------------------------------
    F = __calloc_poly_lst();
    F_ref = __calloc_poly_lst();
    init_poly(f, ctx);
    init_poly(g, ctx);
//-------------------------------------------------------
    #if __DEBUG_F4
        printf("---------------------------------------reduction---------------------------------------\n");
    #endif
    // Формирование "матрицы" F 

    #if __DEBUG_F4
        printf("statr preprocessing\n");
    #endif

    preprocessing(F, Pd, G, ctx);

    #if __DEBUG_F4
        printf("preprocessing completed\n");
    #endif

    
    #if __DEBUG_F4
        printf("start ref\n");
    #endif
    
    // Приведение "матрицы" к верхне треугольному виду 
    ref3(F_ref, F, field, ctx);
    
    #if __DEBUG_F4
        printf("ref completed\n");
    #endif

    
    #if __DEBUG_F4
        printf("select new poly for basis\n");
    #endif

    // while(1){}
    // Выбираем полиномы для добавления в базис

    i = 0;
    while(i < F_ref->len){
        h = g_array_index(F_ref, Polynom, i);
        HM(f, h, ctx);
        flag = 1;

        hg = (Polynom*)F->data;
        for(j = 0; j < F->len; j++){
            HM(g, *hg, ctx);
            if (fq_nmod_mpoly_equal(f, g, ctx) == 1){
                flag = 0;
                break;
            }
            hg++;
        }

        if (flag == 1){
            g_array_append_val(F_, h);
            g_array_remove_index(F_ref, i);
        } else i++;
    }
    #if __DEBUG_F4
        printf("select new poly for basis completed\n");
        printf("---------------------------------------reduction-end---------------------------------------\n");
    #endif
//-------------------------------------------------------
    free_poly_lst(F, ctx);
    free_poly_lst(F_ref, ctx);
    clear_poly(f, ctx);
    clear_poly(g, ctx);
}



// Критерий
void F4_GMI(GArray* P, const GArray* G, const Polynom h, ulong t, const PolynomRing ctx){
    GArray* _P;
    Polynom f, L;
    fq_nmod_mpoly_t hm_h, div, _lcm;
    Polynom* hp;
    F4Pair f4p, f4p_f, f4p_g;
    ulong i, j;
//-------------------------------------------------------
    _P = g_array_new(FALSE, FALSE, sizeof(F4Pair));

    init_poly(div, ctx);
    init_poly(_lcm, ctx);
    init_poly(hm_h, ctx);
    HM(hm_h, h, ctx);

    hp = (Polynom*)G->data;
    for(i = 0; i < t; i++){
        F4Pair new_pair;
        init_F4Pair(&new_pair, *hp, h, ctx);
        g_array_append_val(_P, new_pair);
        hp++;
    }
//-------------------------------------------------------
    // printf("t=%ld\n", t);
    // printf("_P:\n");
    // print_F4Pair_lst(_P, ctx);
    // printf("\n");
    // printf("P:\n");
    // print_F4Pair_lst(P, ctx);
    // printf("\n");

    i = 0;
    while(i < P->len){
        f4p = g_array_index(P, F4Pair, i);

        if (monom_divides(f4p.lcm, hm_h, ctx) == 1){
            lcm(_lcm, h, f4p.f, ctx);
            if (fq_nmod_mpoly_equal(_lcm, f4p.lcm, ctx) == 0){
                lcm(_lcm, h, f4p.g, ctx);
                if (fq_nmod_mpoly_equal(_lcm, f4p.lcm, ctx) == 0){
                    free_F4Pair(&f4p, ctx);
                    g_array_remove_index(P, i);
                    continue;
                }
            }
        }

        i++; 
    }

    i = 0;
    while(i < _P->len){
        f4p_f = g_array_index(_P, F4Pair, i);

        j = 0;
        while(j < _P->len){
            if (i != j){
                f4p_g = g_array_index(_P, F4Pair, j);

                if (monom_divides(f4p_g.lcm, f4p_f.lcm, ctx) == 1){
                    free_F4Pair(&f4p_g, ctx);
                    g_array_remove_index(_P, j);
                    if (j < i) i--;
                    continue;
                }
            }
            j++;
        }
        i++;
    }

    i = 1;
    int flag = 0;
    while(i < _P->len){
        f4p_f = g_array_index(_P, F4Pair, i);

        j = 0;
        flag = 0;
        while(j < i){
            f4p_g = g_array_index(_P, F4Pair, j);
            
            if (fq_nmod_mpoly_equal(f4p_f.lcm, f4p_g.lcm, ctx) == 1){
                flag = 1;
                break;
            }
            j++;
        }

        if (flag){
            g_array_remove_index(_P, i);
            continue;
        }

        i++;
    }

    i = 0;
    while(i < _P->len){
        f4p = g_array_index(_P, F4Pair, i);
        HM(hm_h, f4p.f, ctx);
        HM(_lcm, f4p.g, ctx);
        fq_nmod_mpoly_gcd(div, hm_h, _lcm, ctx);

        if (fq_nmod_mpoly_is_one(div, ctx) == 1){
            free_F4Pair(&f4p, ctx);
            g_array_remove_index(_P, i);
            continue;
        }
        i++;
    }


    // while(_P->len != 0){
    //     f4p = g_array_index(_P, F4Pair, _P->len-1);
    //     g_array_append_val(P, f4p);
    //     g_array_remove_index(_P, _P->len-1);
    // }

    F4Pair* mas = (F4Pair*)_P->data;
    for(i = 0; i < _P->len; i++){
        g_array_append_val(P, mas[i]);
    }

//-------------------------------------------------------
    // g_array_free(_P, TRUE);
    g_array_free(_P, FALSE);
    clear_poly(div, ctx);
    clear_poly(_lcm, ctx);
    clear_poly(hm_h, ctx);
}

// Не работает
void reduced_F4_GMI(GArray* P, GArray* G, GArray* F_, const PolynomRing ctx){
    // Предполается, что P тоже NULL. Первый шаг F4
    if (F_ == NULL){
        reduce_groebner_basis(G, ctx);
        return;
    }

    // printf("F(%d): \n", F_->len);
    // print_poly_lst(F_, ctx);
    // printf("\n\n");

    // Удаляем лишние полиномы из F_
    reduce_groebner_basis(F_, ctx); // Относительно самого себя

    // printf("F(%d) относительно F: \n", F_->len);
    // print_poly_lst(F_, ctx);
    // printf("\n\n");

    // относительно G
    reduce_groebner_basis_relative(F_, G, ctx);

    // printf("F(%d) относительно F и G: \n", F_->len);
    // print_poly_lst(F_, ctx);
    // printf("\n\n");

    // printf("G(%d): \n", G->len);
    // print_poly_lst(G, ctx);
    // printf("\n\n");

    // теперь тоже для G
    int flag;
    Basis Gbasis;
    Basis Fbasis;
    ulong i, j;

    fq_nmod_mpoly_t mi, mj;

    fq_nmod_mpoly_init(mi, ctx);
    fq_nmod_mpoly_init(mj, ctx);


    Fbasis = (Basis)F_->data;
    
    while(!flag){
        Gbasis = (Basis)G->data;
        flag = 1;

        for(i = 0; i < G->len; i++){
            HM(mi, Gbasis[i], ctx);

            // Относительно себя
            for(j = 0; j < G->len; j++){
                if (i == j) continue;
                HM(mj, Gbasis[j], ctx);

                if(monom_divides(mi, mj, ctx)){
                    flag = 0;
                    
                    // Нужно удалить лишние наборы 
                    remove_pairs_containig(P, Gbasis[i], ctx);

                    // Удаляем полином
                    g_array_remove_index(G, i);
                    break;
                }
            }

            if (!flag) break;

            // Относительно F
            for(j = 0; j < F_->len; j++){
                HM(mj, Fbasis[j], ctx);

                if(monom_divides(mi, mj, ctx)){
                    flag = 0;
                    remove_pairs_containig(P, Gbasis[i], ctx);
                    g_array_remove_index(G, i);
                    break;
                }
            }

            if (!flag) break;
        }
    }

    fq_nmod_mpoly_clear(mi, ctx);
    fq_nmod_mpoly_clear(mj, ctx);

    // printf("Конец: ");

    // printf("F(%d): \n", F_->len);
    // print_poly_lst(F_, ctx);
    // printf("\n\n");

    // printf("G(%d): \n", G->len);
    // print_poly_lst(G, ctx);
    // printf("\n\n");

    // Теперь можно использовать GMI и обновить набор пар
    Basis temp = (Basis)F_->data;
    for(i = 0; i < F_->len; i++){
        F4_GMI(P, G, temp[i], G->len, ctx);
        g_array_append_val(G, temp[i]);
    }
}

#ifdef __cplusplus
extern "C"
#endif
F4Result F4(const Basis F, ulong npoly, const Field field, const PolynomRing ctx){
    GArray *G; // Строящийся базис гребнера
    GArray *F_; // Новые полиномы добавляемые в базис
    GArray *P; // Критические пары
    GArray *Pd; // Выбранные критические пары
    ulong d;
    ulong i, j;
    Polynom f, g, h;
    Polynom* hp;
    Basis temp;
//-------------------------------------------------------
    G = __calloc_poly_lst(); // g_array_new(FALSE, FALSE, sizeof(Polynom));
    F_ = __calloc_poly_lst();
    P = g_array_new(FALSE, FALSE, sizeof(F4Pair));
    Pd = g_array_new(FALSE, FALSE, sizeof(F4Pair));

    // printf("LOL");
    fq_nmod_t c;
    fmpz_t t;
    fmpz_t inv_t;
    fmpz_t mod;
    fmpz_init(t);
    fmpz_init(inv_t);
    fmpz_init(mod);
    fq_nmod_init(c, ctx);

    // fmpz_set_ui(mod, fq_nmod_ctx_prime(field));
    fq_nmod_ctx_order(mod, field);
    
    for(i = 0; i < npoly; i++){
        g = __calloc_poly();
        init_poly(g, ctx);
        // set_poly(g, F[i], ctx);

        // Нормируем многочлен g
        HC(c, F[i], ctx);
        fq_nmod_get_fmpz(t, c, ctx);
        fmpz_invmod(inv_t, t, mod);
        fq_nmod_set_fmpz(c, inv_t, ctx);

        // fmpz_print(mod);
        // printf("\n");
        // fmpz_print(t);
        // printf(" inv: ");
        // fmpz_print(inv_t);
        // printf("\n");

        fq_nmod_mpoly_scalar_mul_n_fq(g, F[i], inv_t, ctx);
        g_array_append_val(G, g);
    }

    fmpz_clear(t);
    fmpz_clear(inv_t);
    fmpz_clear(mod);
    fq_nmod_clear(c, ctx);

    // Формирование пар с учетом lcm и gcd критерия
    for(i = 0; i < npoly; i++)
        F4_GMI(P, G, g_array_index(G, Polynom, i), i, ctx);
//-------------------------------------------------------

    #if __DEBUG_F4
        printf("%d\n", G->len);
        printf("G:\n");
        print_poly_lst(G, ctx);
        printf("\n");

        printf("P:\n");
        print_F4Pair_lst(P, ctx);
        printf("\n");

        printf("Pd:\n");
        print_F4Pair_lst(Pd, ctx);
        printf("\n");
    #endif

    int counter = 0;

    while(P->len > 0){
        #if __DEBUG_F4
            d = find_min_deg_in_F4Pairs(P, ctx);
                // if (counter == 2) break;
                printf("min deg=%ld\n", d);
        #endif

        // Выбираем критические пары, переносим их в Pd и удаляем из P
        F4_select(Pd, P, ctx);

        #if __DEBUG_F4
                printf("Pd:\n");
                print_F4Pair_lst(Pd, ctx);
                printf("\nP:\n");
                print_F4Pair_lst(P, ctx);
                printf("\n");
        #endif

        // Строим новые полиномы по критическим парам в Pd и редуцируем их
        reduction(F_, Pd, G, field, ctx);

        #if __DEBUG_F4
            printf("F+:\n");
            print_poly_lst(F_, ctx);
            printf("\n");
        #endif

        
        // reduced_F4_GMI(P, G, F_, ctx);
        temp = (Basis)F_->data;
        for(i = 0; i < F_->len; i++){
            F4_GMI(P, G, temp[i], G->len, ctx);
            g_array_append_val(G, temp[i]);
        }

        g_array_free(F_, FALSE);
        F_ = __calloc_poly_lst();

        // while(F_->len > 0){
        //     f = g_array_index(F_, Polynom, F_->len-1);
            
        //     // Добавляем новые критические пары
        //     F4_GMI(P, G, f, G->len, ctx);

        //     // Добавляем новый полином в базис
        //     g_array_append_val(G, f);
        //     g_array_remove_index(F_, F_->len-1);
        // }

        #if __DEBUG_F4
            printf("G len: %d\n", G->len);
            printf("P len: %d\n", P->len);
            printf("Pd len: %d\n", Pd->len);
        

            printf("P:\n");
            print_F4Pair_lst(P, ctx);
            printf("\n");

            printf("Pd:\n");
            print_F4Pair_lst(Pd, ctx);
            printf("\n");

            printf("G:\n");
            print_poly_lst(G, ctx);
            printf("\n");

            sleep(3);
        #endif
        // counter++;
        // break;
    }

    // #if __DEBUG_CHECK
    //     Basis b = (Basis)G->data;
    //     fq_nmod_mpoly_t s, s_mod;
    //     Basis Q = init_empty_basis(G->len, ctx);
    //     int flag = 1;

    //     init_poly(s, ctx);
    //     init_poly(s_mod, ctx);

    //     for(i = 0; i < G->len-1; i++){
    //         for(j = i+1; j < G->len; j++){
    //             // spol(s, b[i], b[j], field, ctx);
    //             spol_old(s, b[i], b[j], ctx);
    //             fq_nmod_mpoly_divrem_ideal(Q, s_mod, s, b, G->len, ctx);
    //             if (fq_nmod_mpoly_is_zero(s_mod, ctx)){
    //                 flag = 0;
    //                 break;
    //             }
    //         }

    //         if (!flag) break;
    //     }

    //     if (flag) printf("This is Groebner basis :)\n");
    //     else printf("This is not Groebner basis :c\n");

    //     free_basis(Q, G->len, ctx);
    //     clear_poly(s, ctx);
    //     clear_poly(s_mod, ctx);
    // #endif
//-------------------------------------------------------
    reduce_groebner_basis(G, ctx);
    Basis res = from_garray(G);
    F4Result resres = {res, G->len};
    // free_poly_lst(G, ctx);
    g_array_free(G, FALSE);
    free_poly_lst(F_, ctx);
    free_F4Pair_lst(P, ctx);
    free_F4Pair_lst(Pd, ctx);
    // g_array_free(Pd, TRUE);

    return resres;
}
