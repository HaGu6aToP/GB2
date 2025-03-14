#include "headers.h"

// First aurgument is the file name. 
// File consider the number of polynomials, field ordering, number of variables and the variables with the polynomials themselves. 
// Each data on a new line.
// Example:
//  2
//  7
//  3
//  x y z
//  x^3*y^2 - x^2*y^3 + x
//  3*x^4*y + y^2
void main(int argc, char** argv){
    if (argc == 0){
        printf("Not enought parameters");
        return;
    }

    char* file_name = argv[1];
    int repeats = parseInt(argv[2]);
    
    FILE* file;
    if ((file = fopen(file_name, "r")) == NULL){
        printf("Failed to open file");
        return;
    }

    ulong npoli; // Polinomials count
    ulong nvars; // Variables count
    char buff[BUFFER_SIZE];
    const char** variables;
    const ordering_t order = ORD_DEGREVLEX; // Ordering
    ulong p; // Field order
    fq_nmod_ctx_t field_ctx; // Field
    fq_nmod_mpoly_ctx_t poly_ring_ctx; // Ring
    Basis basis;
    
    
    fscanf(file, "%ld\n%ld\n%ld\n", &npoli, &p, &nvars);
    fgets(buff, BUFFER_SIZE, file);

    variables = flint_calloc(nvars, sizeof(char*));
    get_variables(variables, nvars, buff);

    // Field init
    fq_nmod_ctx_init_conway(field_ctx, &p, 1, "x");

    // Ring init
    fq_nmod_mpoly_ctx_init(poly_ring_ctx, nvars, order, field_ctx);

    // Reading polinomials
    basis = init_empty_basis(npoli, poly_ring_ctx);
    read_polinomials(basis, npoli, variables, poly_ring_ctx, file);

    fclose(file);

    //------------------------------execution------------------------------
    
    printf("Basis:\n");
    print_basis(basis, npoli, variables, poly_ring_ctx);


    Buchberger_result GBasis = log_buchberger(basis, npoli, poly_ring_ctx);

    print_basis(GBasis.basis, GBasis.len, variables, poly_ring_ctx);

    int check = is_groebner_basis(GBasis.basis, GBasis.len, poly_ring_ctx);
    if (check == 1) printf("This is Groebner basis :)\n");
    else printf("This is not Groebner basis :c\n");
    
    // ----------------------free resources----------------------
    free_basis(basis, npoli, poly_ring_ctx);
    fq_nmod_mpoly_ctx_clear(poly_ring_ctx);
    fq_nmod_ctx_clear(field_ctx);
    free_variables(variables, nvars);
    
}