#include "sparse_matrix.h"

void sparse_matrix_test(){
    sparse_matrix_t sparse_M;
    fq_nmod_ctx_t field;
    ulong p = 7;
    fq_nmod_ctx_init(field, &p, 1, "X");

    sparse_matrix_init(sparse_M, 4, 5, field);

    // printf("lines-main matrix:\n");

    // // 1, 0, 1, 0, 0
    // // 1, 1, 0, 0, 0
    // // 0, 1, 0, 1, 0
    // // 0, 0, 0, 6, 1


    sparse_matrix_add_elem_ui(sparse_M, 0, 2, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 1, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 3, 1);
    sparse_matrix_add_elem_ui(sparse_M, 3, 3, 6);
    sparse_matrix_add_elem_ui(sparse_M, 3, 4, 1);
    sparse_matrix_add_elem_ui(sparse_M, 0, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 1, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 1, 1, 1);

    sparse_matrix_print(sparse_M);
    printf("\n");

    // sparse_matrix_add_elem(sparse_M, 1, 1, 1);
    
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    
    printf("swap columns, (0, 3), (1, 4)\n");
    sparse_matrix_swap_columns(sparse_M, 0, 3);
    sparse_matrix_swap_columns(sparse_M, 1, 4);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");

    printf("swap lines, (0, 3), (0, 2)");
    sparse_matrix_swap_lines(sparse_M, 0, 3);
    sparse_matrix_swap_lines(sparse_M, 0, 2);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");

    sparse_matrix_clear(sparse_M);

    // 5 4 2
    // 0 4 5
    // 6 0 0
    // 0 1 5
    // 2 0 2
    // 0 3 0

    sparse_matrix_init(sparse_M, 6, 3, field);
    sparse_matrix_add_elem_ui(sparse_M, 0, 1, 4);
    sparse_matrix_add_elem_ui(sparse_M, 4, 0, 2);
    sparse_matrix_add_elem_ui(sparse_M, 0, 2, 2);
    sparse_matrix_add_elem_ui(sparse_M, 0, 0, 5);
    sparse_matrix_add_elem_ui(sparse_M, 1, 1, 4);
    sparse_matrix_add_elem_ui(sparse_M, 1, 2, 5);
    sparse_matrix_add_elem_ui(sparse_M, 2, 0, 6);
    sparse_matrix_add_elem_ui(sparse_M, 3, 1, 1);
    sparse_matrix_add_elem_ui(sparse_M, 3, 2, 5);
    sparse_matrix_add_elem_ui(sparse_M, 4, 2, 2);
    sparse_matrix_add_elem_ui(sparse_M, 5, 1, 3);

    sparse_matrix_print(sparse_M);
    printf("\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    printf("swap columns: (0, 1), (1, 2)\n");
    sparse_matrix_swap_columns(sparse_M, 0, 1);
    sparse_matrix_swap_columns(sparse_M, 1, 2);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    
    printf("swap lines: (0, 2), (3, 5)\n");
    sparse_matrix_swap_lines(sparse_M, 0, 2);
    sparse_matrix_swap_lines(sparse_M, 3, 5);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print(sparse_M);

    sparse_matrix_clear(sparse_M);
    fq_nmod_ctx_clear(field);

}