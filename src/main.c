#define _POSIX_C_SOURCE 199309L
#include "buchberger.h"
#include "basis_tools.h"
#include "tools.h"
#include "f4.h"
#include "config.h"

#include <time.h>

#define __DEBUG_CHECK 1

// First aurgument is the file name. 
// Second - amount of repeats
// File consider the number of polynomials, field ordering, number of variables and the variables with the polynomials themselves. 
// Each data on a new line.
// Example:
//  2
//  7
//  3
//  x y z
//  x^3*y^2 - x^2*y^3 + x
//  3*x^4*y + y^2
int main(int argc, char** argv){
    if (argc == 0){
        printf("Not enought parameters\n");
        return 1;
    }

    char* file_name = argv[1];
    int repeats = parseInt(argv[2]);
    
    FILE* file;
    if ((file = fopen(file_name, "r")) == NULL){
        printf("Failed to open file\n");
        return 2;
    }

    ulong npoly; // Polinomials count
    ulong nvars; // Variables count
    char buff[BUFFER_SIZE];
    const char** variables;
    const ordering_t order = ORD_DEGREVLEX; // Ordering
    ulong p; // Field order
    fq_nmod_ctx_t field_ctx; // Field
    fq_nmod_mpoly_ctx_t poly_ring_ctx; // Ring
    Basis basis;
    int check;
    ulong threads_count = 4;
    NO_OF_IRRED = threads_count/2;

    // NO_OF_IRRED = 4;
    // printf("%ld\n\n", NO_OF_IRRED);
    
    
    fscanf(file, "%ld\n%ld\n%ld\n", &npoly, &p, &nvars);
    fgets(buff, BUFFER_SIZE, file);

    // printf("%ld %ld %ld", npoli, p, nvars);

    variables = flint_calloc(nvars, sizeof(char*));
    get_variables(variables, nvars, buff);

    // Field init
    fq_nmod_ctx_init(field_ctx, &p, 1, "x");
    

    // Ring init
    fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, order, field_ctx);

    // Reading polinomials
    basis = init_empty_basis(npoly, poly_ring_ctx);
    read_polinomials(basis, npoly, variables, poly_ring_ctx, file);

    fclose(file);

    // printf("%d", flint_get_num_available_threads());

    //------------------------------execution------------------------------

    // fq_nmod_mpoly_t f, g, s;
    // fq_nmod_mpoly_init(f, poly_ring_ctx);
    // fq_nmod_mpoly_init(g, poly_ring_ctx);
    // fq_nmod_mpoly_init(s, poly_ring_ctx);

    // const char* arr_variables[3] = {"x", "y", "z"};
    // // const char** variables = arr_variables;

    // fq_nmod_mpoly_set_str_pretty(f, "x^3*y^2 - x^2*y^3 + x", arr_variables, poly_ring_ctx);
    // fq_nmod_mpoly_set_str_pretty(g, "3*x^4*y + y^2", arr_variables, poly_ring_ctx);

    // spol(s, f, g, poly_ring_ctx);
    // fq_nmod_mpoly_print_pretty(s, arr_variables, poly_ring_ctx);
    
    printf("Basis:\n");
    print_basis(basis, npoly, variables, poly_ring_ctx);

     F4Result GBasis;

    #if __DEBUG_CHECK
        GBasis = F4(basis, npoly, field_ctx, poly_ring_ctx);
        printf("Groebner basis(F4):\n");

        // min_groebner_basis(GBasis.basis, &GBasis.len, poly_ring_ctx);

        print_basis(GBasis.basis, GBasis.len, variables, poly_ring_ctx);
        printf("Basis len: %ld\n", GBasis.len);

        FILE* f = fopen("F4_res.txt", "w");
        int maxlen = 0;
        fprintf(f, "Basis len: %ld\n", GBasis.len);
        for(ulong i = 0; i < GBasis.len; ++i){
            // fprintf(f, "len %d: ", fq_nmod_mpoly_length(GBasis.basis[i], poly_ring_ctx));
            if (fq_nmod_mpoly_length(GBasis.basis[i], poly_ring_ctx) > maxlen) maxlen = fq_nmod_mpoly_length(GBasis.basis[i], poly_ring_ctx);
            fq_nmod_mpoly_fprint_pretty(f, GBasis.basis[i], NULL, poly_ring_ctx);
            fprintf(f, "\n");
            
        }
        fprintf(f, "Max poly-len: %d", maxlen);
        fclose(f);
    
        check = is_groebner_basis(GBasis.basis, GBasis.len, poly_ring_ctx);
        if (check == 1) printf("This is Groebner basis :)\n");
        else printf("This is not Groebner basis :c\n");
        free_basis(GBasis.basis, GBasis.len, poly_ring_ctx); 
    #endif

    // printf("\n");


    #if __DEBUG_CHECK
        Buchberger_result GBasis2 = buchberger_v2_1(basis, npoly, poly_ring_ctx);
        // Buchberger_result GBasis2 = buchberger(basis, npoly, poly_ring_ctx);
        // Buchberger_result GBasis = log_threaded_buchberger(basis, npoli, threads_count, poly_ring_ctx);
        // Buchberger_result GBasis = threaded_buchberger_v2(basis, npoli, threads_count, poly_ring_ctx);

        printf("Groebner basis(buchberger):\n");
        // min_groebner_basis(GBasis2.basis, &GBasis2.len, poly_ring_ctx);
        print_basis(GBasis2.basis, GBasis2.len, variables, poly_ring_ctx);
        printf("Basis len: %ld\n", GBasis2.len);

        f = fopen("Buchverger_res.txt", "w");
        fprintf(f, "Basis len: %ld\n", GBasis2.len);
        maxlen = 0;
        for(ulong i = 0; i < GBasis2.len; ++i){
            fprintf(f, "len %d: ", fq_nmod_mpoly_length(GBasis2.basis[i], poly_ring_ctx));
            if (fq_nmod_mpoly_length(GBasis2.basis[i], poly_ring_ctx) > maxlen) maxlen = fq_nmod_mpoly_length(GBasis2.basis[i], poly_ring_ctx);
            fq_nmod_mpoly_fprint_pretty(f, GBasis2.basis[i], NULL, poly_ring_ctx);
            fprintf(f, "\n");
            
        }
        fprintf(f, "Max poly-len: %d", maxlen);
        fclose(f);

        check = is_groebner_basis(GBasis2.basis, GBasis2.len, poly_ring_ctx);
        if (check == 1) printf("This is Groebner basis :)\n");
        else printf("This is not Groebner basis :c\n");
        free_basis(GBasis2.basis, GBasis2.len, poly_ring_ctx);
    #endif


    // fq_nmod_mpoly_t polynom1;
    // fq_nmod_mpoly_init(polynom1, poly_ring_ctx);
    // fq_nmod_mpoly_t polynom2;
    // fq_nmod_mpoly_init(polynom2, poly_ring_ctx);
    // fq_nmod_mpoly_t S_poly;
    // fq_nmod_mpoly_init(S_poly, poly_ring_ctx);

    // fq_nmod_mpoly_set_str_pretty(polynom1, "4*x3^2*x4^5 + 5*x3^2*x4^4 + 6*x3^2*x4^3 + 4*x3^2*x4^2 + 4*x3^2*x4 + 5*x3^2", NULL, poly_ring_ctx);
    // fq_nmod_mpoly_set_str_pretty(polynom2, "x3*x4^5 + 3*x3*x4^4 + 5*x3*x4^3 + x3*x4^2 + x3*x4 + 3*x3", NULL, poly_ring_ctx);

    // print_poly("P1:", polynom1, NULL, poly_ring_ctx);
    // print_poly("P2:", polynom2, NULL, poly_ring_ctx);

    // log_S(S_poly, polynom1, polynom2, poly_ring_ctx);
    // print_poly("S:", S_poly, NULL, poly_ring_ctx);

    // fq_nmod_mpoly_clear(S_poly, poly_ring_ctx);
    // fq_nmod_mpoly_clear(polynom1, poly_ring_ctx);
    // fq_nmod_mpoly_clear(polynom2, poly_ring_ctx);

    
    // log_buchberger_v2(basis, npoli, poly_ring_ctx);

    //------------------------testing------------------------
    struct timespec start, end;
    double summ_time = 0;

    
    for (int i = 0; i < repeats; i++){
        clock_gettime(CLOCK_MONOTONIC, &start);
        Buchberger_result GBasis = buchberger_v2_1(basis, npoly, poly_ring_ctx);
        // threaded_buchberger_v2(basis, npoly, threads_count, poly_ring_ctx);
        // GBasis = F4(basis, npoly, field_ctx, poly_ring_ctx);
        clock_gettime(CLOCK_MONOTONIC, &end);
        // printf("%d\n", GBasis.len);
        free_basis(GBasis.basis, GBasis.len, poly_ring_ctx);
        summ_time += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        if (i % 100 == 0) printf("i=%d\n", i);
    }

    printf("Runnig time: %f s\n", summ_time/repeats);

    // timeit_t t;
    // slong summary_time = 0;
    // for (int i = 0; i < repeats; i++){
    //     timeit_start(t);
    //     buchberger(basis, npoli, poly_ring_ctx);
    //     timeit_stop(t);
    //     summary_time += t->cpu;
    // }

    // printf("CPU time: %ld ms\n", summary_time/repeats);
    
    // ----------------------free resources----------------------
    free_basis(basis, npoly, poly_ring_ctx);
    // free_basis(GBasis.basis, GBasis.len, poly_ring_ctx);
    // free_basis(GBasis2.basis, GBasis2.len, poly_ring_ctx);
    fq_nmod_mpoly_ctx_clear(poly_ring_ctx);
    fq_nmod_ctx_clear(field_ctx);
    free_variables(variables, nvars);
    
    return 0;
}
