#include "headers.h"

// void printF(GArray* F, PolynomRing ctx){
//     Polynom h;
//     printf("F: \n");
//     for (int i = 0; i < F->len; i++){
//         h = g_array_index(F, Polynom, i);
//         fq_nmod_mpoly_print_pretty(h, NULL, ctx);
//         printf("\n");
//     }
//     printf("\n");
// }

// void printP(GArray* P){
//     ulong i, j, k;
//     Pair pair;
//     printf("P:\n");
//     for(k = 0; k < P->len; k++){
//         pair = g_array_index(P, Pair, k);
//         i = pair.first;
//         j = pair.second;
//         printf("(%ld, %ld)\n", i, j);
//     }
//     printf("\n");
// }

// void printSPair(SPair spair, PolynomRing ctx){
//     printf("(");
//     fq_nmod_mpoly_print_pretty(spair.poly, NULL, ctx);
//     printf(", %ld, %ld)\n", spair.first, spair.second);

// }

void min_find_test_v1(){
    // ulong nvars = 3;
//     const ordering_t order = ORD_LEX; // Ordering
//     ulong p = 7; // Field order
//     fq_nmod_ctx_t field_ctx; // Field
//     fq_nmod_mpoly_ctx_t poly_ring_ctx; // Ring
//     Basis basis;
    

//     // Field init
//     fq_nmod_ctx_init(field_ctx, &p, 1, "x");
//     // Ring init
//     fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, order, field_ctx);


//     fq_nmod_mpoly_t p1, p2, p3;
//     fq_nmod_mpoly_init(p1, poly_ring_ctx);
//     fq_nmod_mpoly_init(p2, poly_ring_ctx);
//     fq_nmod_mpoly_init(p3, poly_ring_ctx);
//     Polynom pp1 = p1, pp2 = p2, pp3 = p3;
//     Pair pa1 = {0, 1}, pa2 = {0, 2}, pa3 = {1, 2};

//     GArray* P = g_array_new(FALSE, FALSE, sizeof(Pair));
//     g_array_append_val(P, pa1);
//     g_array_append_val(P, pa2);
//     g_array_append_val(P, pa3);

//     GArray* F = g_array_new(FALSE, FALSE, sizeof(fq_nmod_mpoly_struct));
//     g_array_append_val(F, pp1);
//     g_array_append_val(F, pp2);
//     g_array_append_val(F, pp3);


//     fq_nmod_mpoly_set_str_pretty(p1, "x1", NULL, poly_ring_ctx);
//     fq_nmod_mpoly_set_str_pretty(p2, "x2", NULL, poly_ring_ctx);
//     fq_nmod_mpoly_set_str_pretty(p3, "x3", NULL, poly_ring_ctx);
    
//     printF(F, poly_ring_ctx);
//     printP(P);

//     SPair res; 
//     res = find_min_v1(F, P, poly_ring_ctx);
//     printSPair(res, poly_ring_ctx);
//     printf("---------------------------------\n\n");

//     fq_nmod_mpoly_set_str_pretty(p1, "x1*x2 + x1", NULL, poly_ring_ctx);
//     fq_nmod_mpoly_set_str_pretty(p2, "x2^2 + x3^2", NULL, poly_ring_ctx);
//     fq_nmod_mpoly_set_str_pretty(p3, "x1*x3 + x2*x3", NULL, poly_ring_ctx);
//     printF(F, poly_ring_ctx);
//     printP(P);

//     res = find_min_v1(F, P, poly_ring_ctx);
//     printf("res=");
//     printSPair(res, poly_ring_ctx);

//     printf("-----------------------------------------\n");


//     fq_nmod_mpoly_clear(p1, poly_ring_ctx);
//     fq_nmod_mpoly_clear(p2, poly_ring_ctx);
//     fq_nmod_mpoly_clear(p3, poly_ring_ctx);
//     fq_nmod_mpoly_ctx_clear(poly_ring_ctx);
//     fq_nmod_ctx_clear(field_ctx);
//     g_array_free(P, TRUE);
//     g_array_free(F, TRUE);
    
    
// }

// void min_find_test_v2(){
//     ulong nvars = 3;
//     const ordering_t order = ORD_LEX; // Ordering
//     ulong p = 7; // Field order
//     fq_nmod_ctx_t field_ctx; // Field
//     fq_nmod_mpoly_ctx_t poly_ring_ctx; // Ring
//     Basis basis;
    

//     // Field init
//     fq_nmod_ctx_init(field_ctx, &p, 1, "x");
//     // Ring init
//     fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, order, field_ctx);


//     fq_nmod_mpoly_t p1, p2, p3, s01, s02, s12;
//     fq_nmod_mpoly_init(p1, poly_ring_ctx);
//     fq_nmod_mpoly_init(p2, poly_ring_ctx);
//     fq_nmod_mpoly_init(p3, poly_ring_ctx);
//     fq_nmod_mpoly_init(s01, poly_ring_ctx);
//     fq_nmod_mpoly_init(s02, poly_ring_ctx);
//     fq_nmod_mpoly_init(s12, poly_ring_ctx);

//     Polynom pp1 = p1, pp2 = p2, pp3 = p3;
//     Polynom ps01 = s01, ps02 = s02, ps12 = s12;

//     fq_nmod_mpoly_set_str_pretty(p1, "x1*x2 + x1", NULL, poly_ring_ctx);
//     fq_nmod_mpoly_set_str_pretty(p2, "x2^2 + x3^2", NULL, poly_ring_ctx);
//     fq_nmod_mpoly_set_str_pretty(p3, "x1*x3 + x2*x3", NULL, poly_ring_ctx);

//     S(s01, p1, p2, poly_ring_ctx);
//     S(s02, p1, p3, poly_ring_ctx);
//     S(s12, p2, p3, poly_ring_ctx);

//     SPair sp1 = {s01, 0, 1}, sp2 = {s02, 0, 2}, sp3 = {s12, 1, 2};

//     GArray* P = g_array_new(FALSE, FALSE, sizeof(SPair));
//     g_array_append_val(P, sp1);
//     g_array_append_val(P, sp2);
//     g_array_append_val(P, sp3);

//     int res = find_min(P, poly_ring_ctx);
//     printf("res=%d\n", res);
//     printSPair(g_array_index(P, SPair, res), poly_ring_ctx);

//     fq_nmod_mpoly_clear(p1, poly_ring_ctx);
//     fq_nmod_mpoly_clear(p2, poly_ring_ctx);
//     fq_nmod_mpoly_clear(p3, poly_ring_ctx);
//     fq_nmod_mpoly_clear(s01, poly_ring_ctx);
//     fq_nmod_mpoly_clear(s02, poly_ring_ctx);
//     fq_nmod_mpoly_clear(s12, poly_ring_ctx);
//     fq_nmod_mpoly_ctx_clear(poly_ring_ctx);
//     fq_nmod_ctx_clear(field_ctx);
//     g_array_free(P, TRUE);

}