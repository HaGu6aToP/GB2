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
    printf(", deg=%ld", p->deg);
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

    p->deg = deg(p->lcm, ctx);

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
    // F4Pair f4p;

    // while(P->len > 0){
    //     f4p = g_array_index(P, F4Pair, P->len-1);
    //     free_F4Pair(&f4p, ctx);
    //     g_array_remove_index(P, P->len-1);
    // }

    // g_array_free(P, TRUE);

    F4Pair* mas = (F4Pair*)P->data;
    for(ulong i = 0; i < P->len; i++){
        free_F4Pair(&mas[i], ctx);
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
    ulong d = pf4p->deg; // deg(pf4p->lcm, ctx);
    ulong a = 0;
    pf4p++;

    for(ulong i = 1; i < P->len; i++){
        a = pf4p->deg; // deg(pf4p->lcm, ctx);
        if (a < d) d = a;
        pf4p++;
    }

    return d;
}

void preprocessing(GArray* F, GArray* Pd, const GArray* G, const PolynomRing ctx){

}

void ref(GArray* F_ref, const GArray* F, const Field field, const PolynomRing ctx){

}

void reduction(GArray* F_, GArray* Pd, const GArray* G, const Field field, const PolynomRing ctx){

}

void F4_GMI(GArray* P, const GArray* G, const Polynom h, ulong t, const PolynomRing ctx){
    GArray* _P = g_array_new(FALSE, FALSE, sizeof(F4Pair));
    Polynom f, L;
    fq_nmod_mpoly_t hm_h, div, _lcm;
    Polynom* hp;
    F4Pair f4p, f4p_f, f4p_g;
    Basis b = (Basis)G->data;
    F4Pair* mas;
    ulong i, j;

    init_poly(div, ctx);
    init_poly(_lcm, ctx);
    init_poly(hm_h, ctx);
    HM(hm_h, h, ctx);

    for(i = 0; i < t; i++){
        F4Pair pf4p = {};
        init_F4Pair(&pf4p, b[i], h, ctx);
        g_array_append_val(_P, pf4p);
    }

    #if __DEBUG_F4
        printf("t=%ld\n", t);
        printf("_P:\n");
        print_F4Pair_lst(_P, ctx);
        printf("\n");
        printf("P:\n");
        print_F4Pair_lst(P, ctx);
        printf("\n");
    #endif

    // lcm критерий
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

    // lcm критерий
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

    // gcd критерий
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

    mas = (F4Pair*)_P->data;
    for(i = 0; i < _P->len; i++)
        g_array_append_val(P, mas[i]);

    #if __DEBUG_F4
        printf("res P:\n");
        print_F4Pair_lst(P, ctx);
        printf("\n");
    #endif
    

    g_array_free(_P, TRUE);
    clear_poly(div, ctx);
    clear_poly(_lcm, ctx);
    clear_poly(hm_h, ctx);
}

F4Result F4(const Basis F, ulong npoly, const Field field, const PolynomRing ctx){
    GArray *G; // Строящийся базис гребнера
    GArray *F_; // Новые полиномы добавляемые в базис
    GArray *P; // Критические пары
    GArray *Pd; // Выбранные критические пары
    ulong d;
    ulong i, j;
    Polynom f, g, h;
    Polynom* hp;

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

    while(P->len > 0){
        
    }

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
