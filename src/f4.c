#include "f4.h"
#include "tools.h"

#define HM(res, f, ctx) fq_nmod_mpoly_get_term_monomial(res, f, 0, ctx)
#define HT(res, f, ctx) fq_nmod_mpoly_get_term(res, f, 0, ctx)
#define HC(res, f, ctx) fq_nmod_mpoly_get_term_coeff_fq_nmod(res, f, 0, ctx)

#define init_poly(f, ctx) fq_nmod_mpoly_init(f, ctx)
#define clear_poly(f, ctx) fq_nmod_mpoly_clear(f, ctx)
#define set_poly(res, f, ctx) fq_nmod_mpoly_set(res, f, ctx)

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

    // ulong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    // fq_nmod_mpoly_t lcm_poly;
    // fq_nmod_mpoly_init(lcm_poly, ctx);

    // fq_nmod_mpoly_one(lcm_poly, ctx);
    // lcm(lcm_poly, p1, p2, ctx);

    // fq_nmod_mpoly_t leading_monom_p1, leading_monom_p2, A;
    // fq_nmod_mpoly_init(leading_monom_p1, ctx);
    // fq_nmod_mpoly_init(leading_monom_p2, ctx);
    // fq_nmod_mpoly_init(A, ctx);

    // fq_nmod_mpoly_get_term(leading_monom_p1, p1, 0, ctx);
    // fq_nmod_mpoly_get_term(leading_monom_p2, p2, 0, ctx);

    // fq_nmod_mpoly_div(A, lcm_poly, leading_monom_p1, ctx); // LCM(p1, p2) / LT(f)
    // fq_nmod_mpoly_mul(leading_monom_p1, A, p1, ctx); // LCM(p1, p2) / LT(f) * p1

    // fq_nmod_mpoly_div(A, lcm_poly, leading_monom_p2, ctx); //LCM(p1, p2) / LT(f)
    // fq_nmod_mpoly_mul(leading_monom_p2, A, p2, ctx); //LCM(p1, p2) / LT(f) * p2

    // fq_nmod_mpoly_sub(res, leading_monom_p1, leading_monom_p2, ctx); // S

    // fq_nmod_mpoly_clear(lcm_poly, ctx);
    // fq_nmod_mpoly_clear(leading_monom_p1, ctx);
    // fq_nmod_mpoly_clear(leading_monom_p2, ctx);
    // fq_nmod_mpoly_clear(A, ctx);
}