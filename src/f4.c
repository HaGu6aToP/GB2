#include "f4.h"
#include "tools.h"
#include "basis_tools.h"
#include <unistd.h>
#include "flint/fq_nmod_mat.h"

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
    printf(" }\n");
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
}

void print_F4Pair_lst(const GArray* P, const PolynomRing ctx){
    for(ulong i = 0; i < P->len; i++){
        print_F4Pair(&g_array_index(P, F4Pair, i), ctx);
    }
    printf("\n");
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
    while(Pd->len > 0){
        f4p = g_array_index(Pd, F4Pair, Pd->len-1);
        new_poly = __calloc_poly();
        fq_nmod_mpoly_mul(new_poly, f4p.t_f, f4p.f, ctx);
        g_array_append_val(F, new_poly);

        new_poly = __calloc_poly();
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
    printf("---------------------------------------preprocessing---------------------------------------\n");
    
    while(sub->len != 0){
        printf("F:\n");
        print_poly_lst(F, ctx);
        printf("\n");

        printf("Done:\n");
        print_poly_lst(done, ctx);
        printf("\n");

        printf("HM(F)\\Done:\n");
        print_poly_lst(sub, ctx);
        printf("\n");

        k = max_poly_in_lst(sub, ctx);
        printf("k=%ld\n", k);

        f = g_array_index(sub, Polynom, k);
        printf("selected monom: ");
        fq_nmod_mpoly_print_pretty(f, NULL, ctx);
        printf("\n");

        g_array_append_val(done, f);
        g_array_remove_index(sub, k);

        hp = (Polynom*)G->data;
        for(i = 0; i < G->len; i++){
            HM(m, *hp, ctx);
            if (fq_nmod_mpoly_divides(div, f, m, ctx) == 1){
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

                new_poly = __calloc_poly();
                init_poly(new_poly, ctx);
                fq_nmod_mpoly_mul(new_poly, div, *hp, ctx);
                g_array_append_val(F, new_poly);
                printf("new poly - ");
                fq_nmod_mpoly_print_pretty(new_poly, NULL, ctx);
                printf("\n");

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

    printf("F:\n");
    print_poly_lst(F, ctx);
    printf("\n");

    printf("Done:\n");
    print_poly_lst(done, ctx);
    printf("\n");
    printf("HM(F)\\Done:\n");
    print_poly_lst(sub, ctx);
    printf("\n");

    printf("---------------------------------------preprocessing-end---------------------------------------\n");
//-------------------------------------------------------
    free_poly_lst(done, ctx);
    free_poly_lst(sub, ctx);
    clear_poly(m, ctx);
    clear_poly(div, ctx);
}

void ref(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){
    GArray* F_monoms;
    Polynom f;
    Polynom *hp, *hg;
    fq_nmod_mpoly_t m, sum, g;
    fq_nmod_mat_t M;
    fq_nmod_t x;
    ulong i, j, k;
//-------------------------------------------------------
    printf("---------------------------------------ref---------------------------------------\n");
    F_monoms = __calloc_poly_lst();
    monom_lst_from_poly_lst(F_monoms, F, ctx);
    init_poly(m, ctx);
    init_poly(sum, ctx);
    init_poly(g, ctx);


    poly_quick_sort(F_monoms, 0, F_monoms->len-1, 1, ctx);

    fq_nmod_mat_init(M, F->len, F_monoms->len, field);
    fq_nmod_init(x, field);

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
//-------------------------------------------------------
    printf("F:\n");
    print_poly_lst(F, ctx);
    printf("\n");
    printf("F_monoms:\n");
    print_poly_lst(F_monoms, ctx);
    printf("\n");
    printf("M:\n");
    fq_nmod_mat_print_pretty(M, field);
    printf("\n");

    slong* p = flint_calloc(F_monoms->len, sizeof(slong));
    for(i = 0; i < F_monoms->len; i++)
        p[i] = i;
    slong r;

    r = fq_nmod_mat_lu(p, M, 0, field);
    printf("new columns ordering: \n");
    for(i = 0; i < F_monoms->len; i++)
        printf("%ld ", p[i]);
    printf("\n");

    printf("\n");
    printf("M LU:\n");
    fq_nmod_mat_print_pretty(M, field);
    printf("\n");
    printf("rank=%ld\n", r);

    fq_nmod_t y;
    fq_nmod_init(y, field);

    // for(i = 0; i < r; i++){
    //     if (fq_nmod_is_one(fq_nmod_mat_entry(M, i, i), field) == 0) {

    //         fq_nmod_inv(y, fq_nmod_mat_entry(M, i, i), field);
    //         for(j; j < F_monoms->len; j++){
    //             if (fq_nmod_is_zero((const fq_nmod_struct*)fq_nmod_mat_entry(M, i, j), field) == 1) continue;
    //             fq_nmod_mul(x ,y, fq_nmod_mat_entry(M, i, j), field);
    //             fq_nmod_mat_entry_set(M, i, j, x, field);
    //         }
    //     }
    // }

    printf("\n");
    printf("M LU:\n");
    fq_nmod_mat_print_pretty(M, field);
    printf("\n");


    for(i = 0; i < r; i++){
        f = __calloc_poly();
        init_poly(f, ctx);
        
        fq_nmod_mpoly_scalar_mul_fq_nmod(f, g_array_index(F_monoms, Polynom, p[i]), (const fq_nmod_struct*)fq_nmod_mat_entry(M, i, i), ctx);
        fq_nmod_mpoly_zero(sum, ctx);

        for(j = i+1; j < F_monoms->len; j++){
            if (fq_nmod_is_zero((const fq_nmod_struct*)fq_nmod_mat_entry(M, i, j), field) == 1) continue;
            fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, p[j]), (const fq_nmod_struct*)fq_nmod_mat_entry(M, i, j), ctx);
            fq_nmod_mpoly_add(sum, f, m, ctx);
            set_poly(f, sum, ctx);
        }

        g_array_append_val(F_ref, f);
    }

    printf("---------------------------------------ref-end---------------------------------------\n");
//-------------------------------------------------------
    fq_nmod_clear(y, field);
    free_poly_lst(F_monoms, ctx);
    clear_poly(m, ctx);
    clear_poly(sum, ctx);
    clear_poly(g, ctx);
    flint_free(p);
    fq_nmod_mat_clear(M, field);
    fq_nmod_clear(x, field);

}

void reduction(GArray* F_, GArray* Pd, const GArray* G, const Field field, const PolynomRing ctx){
    GArray* F;
    GArray* F_hm;
    GArray* F_ref;
    Polynom *hg;
    Polynom h;
    fq_nmod_mpoly_t f, g;
    ulong i, j;
    int flag;
//-------------------------------------------------------
    // F = g_array_new(FALSE, FALSE, sizeof(Polynom));
    F =__calloc_poly_lst();
    F_ref = __calloc_poly_lst();
    init_poly(f, ctx);
    init_poly(g, ctx);
//-------------------------------------------------------
    printf("---------------------------------------reduction---------------------------------------\n");
    preprocessing(F, Pd, G, ctx);
    ref(F_ref, F, field, ctx);

    printf("F:\n");
    print_poly_lst(F, ctx);
    printf("\nF_ref:\n");
    print_poly_lst(F_ref, ctx);
    printf("\n");

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

    printf("---------------------------------------reduction-end---------------------------------------\n");
//-------------------------------------------------------
    free_poly_lst(F, ctx);
    free_poly_lst(F_ref, ctx);
    clear_poly(f, ctx);
    clear_poly(g, ctx);
}


F4Result F4(const Basis F, ulong npoly, const Field field, const PolynomRing ctx){
    GArray *G, *F_, *P, *Ld, *Pd;
    ulong d;
    ulong i, j;
    Polynom f, g, h;
    Polynom* hp;
//-------------------------------------------------------
    G = g_array_new(FALSE, FALSE, sizeof(Polynom));
    F_ = g_array_new(FALSE, FALSE, sizeof(Polynom));
    P = g_array_new(FALSE, FALSE, sizeof(F4Pair));
    Pd = g_array_new(FALSE, FALSE, sizeof(F4Pair));
    // Ld = g_array_new(FALSE, FALSE, sizeof(F4PairProjection));

    for(i = 0; i < npoly; i++){
        f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
        init_poly(f, ctx);
        init_poly(g, ctx);
        set_poly(f, F[i], ctx);
        set_poly(g, F[i], ctx);

        g_array_append_val(G, g);
        // g_array_append_val(F_, f);
    }

    for(i = 0; i < npoly - 1; i++)
        for(j = i + 1; j < npoly; j++){
            F4Pair p;
            init_F4Pair(&p, F[i], F[j], ctx);
            g_array_append_val(P, p);
        }
//-------------------------------------------------------
        while(P->len > 0){
            printf("G:\n");
            print_poly_lst(G, ctx);
            printf("\n");

            printf("P:\n");
            print_F4Pair_lst(P, ctx);
            printf("\n");

            d = find_min_deg_in_F4Pairs(P, ctx);
            printf("min deg=%ld\n", d);

            free_F4Pair_lst(Pd, ctx);
            F4_select(Pd, P, ctx);

            printf("Pd:\n");
            print_F4Pair_lst(Pd, ctx);
            printf("P:\n");
            print_F4Pair_lst(P, ctx);

            reduction(F_, Pd, G, field, ctx);

            printf("F+:\n");
            print_poly_lst(F_, ctx);
            printf("\n");

            while(F_->len > 0){
                f = g_array_index(F_, Polynom, F_->len-1);
                hp = (Polynom*)G->data;
                for(i = 0; i < G->len; i++){
                    F4Pair new_pair;
                    g = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
                    h = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
                    init_poly(g, ctx);
                    init_poly(h, ctx);
    
                    set_poly(h, f, ctx);
                    set_poly(g, *hp, ctx);
    
                    init_F4Pair(&new_pair, h, g, ctx);
    
                    hp++;
                }
                g_array_append_val(G, h);
                g_array_remove_index(F_, F_->len-1);
            }

            // break;

        }
//-------------------------------------------------------
    Basis res = from_garray(G);
    F4Result resres = {res, G->len};
    return resres;
    }