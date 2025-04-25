#include "f4.h"
#include <unistd.h>
#include "flint/fq_nmod_mat.h"

#include "tools.h"
#include "basis_tools.h"
#include "sparse_matrix.h"

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





void* __calloc_poly_lst(){
    return g_array_new(FALSE, FALSE, sizeof(Polynom));
}
void* __calloc_poly(){
    return flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
}

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

void preprocessing(GArray* F, GArray* Pd, const GArray* G, const PolynomRing ctx){
    F4Pair f4p;
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

    // Формируем *S-пары*
    while(Pd->len > 0){
        f4p = g_array_index(Pd, F4Pair, Pd->len-1);
        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, f4p.t_f, f4p.f, ctx);
        g_array_append_val(F, new_poly);

        new_poly = __calloc_poly();
        init_poly(new_poly, ctx);
        fq_nmod_mpoly_mul(new_poly, f4p.t_g, f4p.g, ctx);
        g_array_append_val(F, new_poly);

        free_F4Pair(&f4p, ctx);
        g_array_pop(Pd);
    }

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
    // printf("---------------------------------------preprocessing---------------------------------------\n");
    
    // Добавляем новые полиномы до тех пор, пока для каждого монома из T(F) не найдется полином f из F для котрого он ведущий
    while(sub->len != 0){
        // printf("F:\n");
        // print_poly_lst(F, ctx);
        // printf("\n");

        // printf("Done:\n");
        // print_poly_lst(done, ctx);
        // printf("\n");

        // printf("HM(F)\\Done:\n");
        // print_poly_lst(sub, ctx);
        // printf("\n");

        k = max_poly_in_lst(sub, ctx);
        // printf("k=%ld\n", k);

        f = g_array_index(sub, Polynom, k);
        // printf("selected monom: ");
        // fq_nmod_mpoly_print_pretty(f, NULL, ctx);
        // printf("\n");

        g_array_append_val(done, f);
        g_array_remove_index(sub, k);
        
        hp = (Polynom*)G->data;
        for(i = 0; i < G->len; i++){
            HM(m, *hp, ctx);
            if (fq_nmod_mpoly_divides(div, f, m, ctx) == 1){
                // printf("selected monom - ");
                // fq_nmod_mpoly_print_pretty(f, NULL, ctx);
                // printf(" divides by HT(");
                // fq_nmod_mpoly_print_pretty(*hp, NULL, ctx);
                // printf(")=");
                // fq_nmod_mpoly_print_pretty(m, NULL, ctx);
                // printf("\n");
                // printf("div=");
                // fq_nmod_mpoly_print_pretty(div, NULL, ctx);
                // printf("\n");

                new_poly = __calloc_poly();
                init_poly(new_poly, ctx);
                fq_nmod_mpoly_mul(new_poly, div, *hp, ctx);
                g_array_append_val(F, new_poly);

                // printf("new poly - ");
                // fq_nmod_mpoly_print_pretty(new_poly, NULL, ctx);
                // printf("\n");

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
    // printf("---------------------------------------preprocessing---------------------------------------\n");
//-------------------------------------------------------
    // g_array_free(sub, TRUE);
    free_poly_lst(sub, ctx);
    free_poly_lst(done, ctx);
    clear_poly(m, ctx);
    clear_poly(div, ctx);
}

// void ref(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){
//     GArray* F_monoms;
//     Polynom f;
//     Polynom *hp, *hg;
//     fq_nmod_mpoly_t m, sum, g;
//     fq_nmod_mat_t M;
//     fq_nmod_t x, y;
//     ulong i, j, k;
//     slong* p;
//     slong r, p_len, t;
// //-------------------------------------------------------
//     F_monoms = __calloc_poly_lst();
//     monom_lst_from_poly_lst(F_monoms, F, ctx);
//     init_poly(m, ctx);
//     init_poly(sum, ctx);
//     init_poly(g, ctx);
//     p_len = MAX(F_monoms->len, F->len);
//     p = flint_calloc(p_len, sizeof(slong));
//     // for(i = 0; i < p_len; i++)
//     //     p[i] = i;

//     poly_quick_sort(F_monoms, 0, F_monoms->len-1, 1, ctx);

//     // printf("%d %d", F->len, F_monoms->len);

//     fq_nmod_mat_init(M, F->len, F_monoms->len, field);
//     fq_nmod_init(x, field);
//     fq_nmod_init(y, field);

//     //     // Формируем матрицу 
//     // hp = (Polynom*)F->data;
//     // for(i = 0; i < F->len; i++){
//     //     for(j = 0; j < fq_nmod_mpoly_length(*hp, ctx); j++){
//     //         fq_nmod_mpoly_get_term_monomial(m, *hp, j, ctx);
//     //         fq_nmod_mpoly_get_term_coeff_fq_nmod(x, *hp, j, ctx);
            
//     //         hg = (Polynom*)F_monoms->data;
//     //         for(k = 0; k < F_monoms->len; k++){
//     //             if (fq_nmod_mpoly_equal(m, *hg, ctx) == 1){
//     //                 fq_nmod_mat_entry_set(M, i, k, x, field);
//     //                 break;
//     //             }
//     //             hg++;
//     //         }
//     //     }
//     //     hp++;
//     // }

//     // printf("M:\n");
//     // fq_nmod_mat_print_pretty(M, field);
//     // printf("\n");

//     // r = fq_nmod_mat_lu_classical(p, M, 0, field);

//     // printf("%d %d\n", F->len, F_monoms->len);
//     // printf("M LU rank=%ld:\n", r);
//     // fq_nmod_mat_print_pretty(M, field);
//     // printf("\n");
//     // // flint_free(p);

//     // slong tt;



//     //     // Получаем редуцированные полиномы
//     // for(i = 0; i < r; i++){
//     //     f = __calloc_poly();
//     //     init_poly(f, ctx);

//     //     for(j=i; j < F_monoms->len; j++){
//     //         if (fq_nmod_is_zero(fq_nmod_mat_entry(M, i, j), field) == 1) continue;

//     //         if (F_monoms->len >= F->len && p[p_len - 1] != -1) tt = p[j]; 
//     //         else tt = j;

//     //         // printf("k=%ld\n", tt);

//     //         fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, tt), fq_nmod_mat_entry(M, i, j), ctx);
//     //         set_poly(sum, f, ctx);
//     //         fq_nmod_mpoly_add(f, sum, m, ctx);
//     //     }

//     //     // g_array_append_val(F_ref, f);
//     //     fq_nmod_mpoly_print_pretty(f, NULL, ctx);
//     //     printf("\n");
//     //     // clear_poly(f, ctx);
//     //     // flint_free(f);
//     // }

// //-------------------------------------------------------

//     // printf("F:\n");
//     // print_poly_lst(F, ctx);
//     // printf("\n");
//     // printf("F_monoms:\n");
//     // print_poly_lst(F_monoms, ctx);
//     // printf("\n");

//     sparse_matrix_t sparse_M;
//     sparse_matrix_init(sparse_M, 4, 5, field);

    

//     // sparse_matrix_add_elem(sparse_M, 0, 2, 1);
//     // sparse_matrix_add_elem(sparse_M, 2, 1, 1);
//     // sparse_matrix_add_elem(sparse_M, 2, 3, 1);
//     // sparse_matrix_add_elem(sparse_M, 3, 3, 6);
//     // sparse_matrix_add_elem(sparse_M, 3, 4, 1);
//     // sparse_matrix_add_elem(sparse_M, 0, 0, 1);
//     // sparse_matrix_add_elem(sparse_M, 1, 0, 1);
//     // sparse_matrix_add_elem(sparse_M, 1, 1, 1);
    
//     // sparse_matrix_print_info(sparse_M);
//     // printf("\n\n");
//     // sparse_matrix_print_pretty(sparse_M);
//     // printf("\n");

//     // sparse_matrix_print(sparse_M);
//     // printf("\n");
//     // sparse_matrix_swap_columns(sparse_M, 0, 3);
//     // sparse_matrix_swap_columns(sparse_M, 1, 4);
//     // // sparse_matrix_canonize(sparse_M);
//     // sparse_matrix_print(sparse_M);
//     // printf("\n");

//     // sparse_matrix_print_pretty(sparse_M);
//     // printf("\n");
//     // sparse_matrix_print_info(sparse_M);
//     // printf("\n\n");
    


//     // sparse_matrix_rem_item(sparse_M, 1, 1);
 
//     // sparse_matrix_print(sparse_M);
//     // printf("\n");
    
    
// //-------------------------------------------------------
//     free_poly_lst(F_monoms, ctx);
//     clear_poly(m, ctx);
//     clear_poly(sum, ctx);
//     clear_poly(g, ctx);
//     fq_nmod_clear(x, field);
//     fq_nmod_clear(y, field);
//     fq_nmod_mat_clear(M, field);
//     flint_free(p);
//     sparse_matrix_clear(sparse_M);
// }


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
    // printf("---------------------------------------ref---------------------------------------\n");
    // printf("F:\n");
    // print_poly_lst(F, ctx);
    // printf("\n");
    // printf("F_monoms:\n");
    // print_poly_lst(F_monoms, ctx);
    // printf("\n");
    // printf("%d %d\n", F_monoms->len, fq_nmod_mat_ncols(M, field));
    // printf("Columns=%d, Lines=%d\n", F_monoms->len, F->len);
    // printf("M:\n");
    // fq_nmod_mat_print_pretty(M, field);
    // printf("\n");

    // printf("columns ordering: \n");
    // for(i = 0; i < p_len; i++)
    //     printf("%ld ", p[i]);
    // printf("\n");

    // Приводим к верхне треугольней форме с помощью LU разложения
    p[p_len - 1] = -1;
    r = fq_nmod_mat_lu_classical(p, M, 0, field);

    // printf("new columns ordering: \n");
    // for(i = 0; i < p_len; i++)
    //     printf("%ld ", p[i]);
    // printf("\n");

    // printf("\n");
    // printf("%d %d\n", F_monoms->len, fq_nmod_mat_ncols(M, field));
    // printf("M LU:\n");
    // fq_nmod_mat_print_pretty(M, field);
    // printf("\n");
    // printf("rank=%ld\n", r);

    // printf("\n");
    // printf("M LU:\n");
    // fq_nmod_mat_print_pretty(M, field);
    // printf("\n");

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
    }

    // printf("---------------------------------------ref-end---------------------------------------\n");
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
    // printf("---------------------------------------reduction---------------------------------------\n");
    // Формирование "матрицы" F 
    preprocessing(F, Pd, G, ctx);

    // Приведение "матрицы" к верхне треугольному виду 
    ref(F_ref, F, field, ctx);


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
    // printf("---------------------------------------reduction-end---------------------------------------\n");
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

        if (fq_nmod_mpoly_divides(div, f4p.lcm, hm_h, ctx) == 1){
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

                if (fq_nmod_mpoly_divides(div, f4p_g.lcm, f4p_f.lcm, ctx) == 1){
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


    while(_P->len != 0){
        f4p = g_array_index(_P, F4Pair, _P->len-1);
        g_array_append_val(P, f4p);
        g_array_remove_index(_P, _P->len-1);
    }

//-------------------------------------------------------
    g_array_free(_P, TRUE);
    clear_poly(div, ctx);
    clear_poly(_lcm, ctx);
    clear_poly(hm_h, ctx);
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
//-------------------------------------------------------
    G = __calloc_poly_lst(); // g_array_new(FALSE, FALSE, sizeof(Polynom));
    F_ = __calloc_poly_lst();
    P = g_array_new(FALSE, FALSE, sizeof(F4Pair));
    Pd = g_array_new(FALSE, FALSE, sizeof(F4Pair));

    for(i = 0; i < npoly; i++){
        g = __calloc_poly();
        init_poly(g, ctx);
        set_poly(g, F[i], ctx);
        g_array_append_val(G, g);
    }

    // Формирование пар с учетом lcm и gcd критерия
    for(i = 0; i < npoly; i++)
        F4_GMI(P, G, g_array_index(G, Polynom, i), i, ctx);
//-------------------------------------------------------

    // printf("%d\n", G->len);
    // printf("G:\n");
    // print_poly_lst(G, ctx);
    // printf("\n");

    // printf("P:\n");
    // print_F4Pair_lst(P, ctx);
    // printf("\n");

    // printf("Pd:\n");
    // print_F4Pair_lst(Pd, ctx);
    // printf("\n");
    // int counter = 0;

    while(P->len > 0){
        // if (counter == 4) break;

        // d = find_min_deg_in_F4Pairs(P, ctx);
        // printf("min deg=%ld\n", d);

        // Выбираем критические пары, переносим их в Pd и удаляем из P
        F4_select(Pd, P, ctx);

        // printf("Pd:\n");
        // print_F4Pair_lst(Pd, ctx);
        // printf("\nP:\n");
        // print_F4Pair_lst(P, ctx);
        // printf("\n");

        // Строим новые полиномы по критическим парам в Pd и редуцируем их
        reduction(F_, Pd, G, field, ctx);

        // printf("F+:\n");
        // print_poly_lst(F_, ctx);
        // printf("\n");

        
        while(F_->len > 0){
            f = g_array_index(F_, Polynom, F_->len-1);
            
            // Добавляем новые критические пары
            F4_GMI(P, G, f, G->len, ctx);

            // Добавляем новый полином в базис
            g_array_append_val(G, f);
            g_array_remove_index(F_, F_->len-1);
        }

        // printf("G len: %d\n", G->len);
        // printf("P len: %d\n", P->len);
        // printf("Pd len: %d\n", Pd->len);

        // printf("P:\n");
        // print_F4Pair_lst(P, ctx);
        // printf("\n");

        // printf("Pd:\n");
        // print_F4Pair_lst(Pd, ctx);
        // printf("\n");

        // printf("G:\n");
        // print_poly_lst(G, ctx);
        // printf("\n");

        // sleep(3);
        // counter++;
        // break;
    }
//-------------------------------------------------------
    Basis res = from_garray(G);
    F4Result resres = {res, G->len};
    // free_poly_lst(G, ctx);
    g_array_free(G, TRUE);
    free_poly_lst(F_, ctx);
    free_F4Pair_lst(P, ctx);
    free_F4Pair_lst(Pd, ctx);
    // g_array_free(Pd, TRUE);

    return resres;
}