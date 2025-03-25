#include "headers.h"

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
