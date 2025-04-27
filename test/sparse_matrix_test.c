#include "sparse_matrix.h"

void sparse_matrix_test(){
    sparse_matrix_t sparse_M;
    fq_nmod_ctx_t field;
    ulong p = 7;
    ulong r;
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
    sparse_matrix_print(sparse_M);
    printf("\n");
    
    printf("swap columns, (0, 3), (1, 4)\n");
    sparse_matrix_swap_columns(sparse_M, 0, 3);
    sparse_matrix_swap_columns(sparse_M, 1, 4);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");
 

    printf("swap lines, (0, 3), (0, 2)\n");
    sparse_matrix_swap_lines(sparse_M, 0, 3);
    sparse_matrix_swap_lines(sparse_M, 0, 2);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("add to line 0 line 1\n");
    sparse_matrix_add_line_mul_ui(sparse_M, 0, 1, 1);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("add to line 1 line 0\n");
    sparse_matrix_add_line_mul_ui(sparse_M, 1, 0, 1);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    printf("add to line 0 line 3 mul 6\n");
    sparse_matrix_add_line_mul_ui(sparse_M, 0, 3, 6);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("mul line 0 by 0\n");
    sparse_matrix_mul_line_ui(sparse_M, 0, 0);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("mul line 1 by 7\n");
    sparse_matrix_mul_line_ui(sparse_M, 1, 7);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");


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

    printf("swap columns: (0, 1)\n");
    sparse_matrix_swap_columns(sparse_M, 0, 1);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");
    

    printf("swap columns: (1, 2)\n");
    sparse_matrix_swap_columns(sparse_M, 1, 2);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");
    
    
    printf("swap lines: (0, 2), (3, 5)\n");
    sparse_matrix_swap_lines(sparse_M, 0, 2);
    sparse_matrix_swap_lines(sparse_M, 3, 5);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("add to line 0 line 1\n");
    sparse_matrix_add_line_mul_ui(sparse_M, 0, 1, 1);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("add to line 3 line 4\n");
    sparse_matrix_add_line_mul_ui(sparse_M, 3, 4, 1);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("add to line 5 line 3 mul 2\n");
    sparse_matrix_add_line_mul_ui(sparse_M, 5, 3, 2);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    printf("mul line 0 by 0\n");
    sparse_matrix_mul_line_ui(sparse_M, 0, 0);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    printf("mul line 1 by 7\n");
    sparse_matrix_mul_line_ui(sparse_M, 1, 7);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    // row echelon form
    // 1 0 0
    // 0 1 0
    // 0 0 1
    printf("gauss ref:\n");
    sparse_matrix_gauss_ref(sparse_M);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    // sparse_matrix_add_line_mul_ui(sparse_M, 3, 0, 5);

    // sparse_matrix_print_info(sparse_M);
    // printf("\n\n");
    // sparse_matrix_print_pretty(sparse_M);
    // printf("\n");
    // sparse_matrix_print_info(sparse_M);
    // printf("\n\n");
    // sparse_matrix_print(sparse_M);
    // printf("\n");

    sparse_matrix_clear(sparse_M);

// 1, 0, 0, 0, 0, 0, 0, 1, 1, 4, 0, 0
// 1, 6, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0
// 0, 1, 6, 1, 1, 1, 4, 0, 0, 0, 0, 0
// 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 6, 2

    sparse_matrix_init(sparse_M, 4, 12, field);

    sparse_matrix_add_elem_ui(sparse_M, 0, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 0, 7, 1);
    sparse_matrix_add_elem_ui(sparse_M, 0, 8, 1);
    sparse_matrix_add_elem_ui(sparse_M, 0, 9, 4);
    
    sparse_matrix_add_elem_ui(sparse_M, 1, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 1, 1, 6);
    sparse_matrix_add_elem_ui(sparse_M, 1, 3, 2);
    
    sparse_matrix_add_elem_ui(sparse_M, 2, 1, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 2, 6);
    sparse_matrix_add_elem_ui(sparse_M, 2, 3, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 4, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 5, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 6, 4);

    sparse_matrix_add_elem_ui(sparse_M, 3, 5, 1);
    sparse_matrix_add_elem_ui(sparse_M, 3, 10, 6);
    sparse_matrix_add_elem_ui(sparse_M, 3, 11, 2);

    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    printf("gauss ref:\n");
    sparse_matrix_gauss_ref(sparse_M);

    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");
    sparse_matrix_print_info(sparse_M);
    printf("\n\n");
    sparse_matrix_print(sparse_M);
    printf("\n");

    sparse_matrix_clear(sparse_M);


    // 1 0 1 0 0 
    // 1 1 0 0 0 
    // 0 1 0 1 0 
    // 0 0 0 6 1

    sparse_matrix_init(sparse_M, 4, 5, field);
    sparse_matrix_add_elem_ui(sparse_M, 0, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 0, 2, 1);

    sparse_matrix_add_elem_ui(sparse_M, 1, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 1, 1, 1);

    sparse_matrix_add_elem_ui(sparse_M, 2, 1, 1);
    sparse_matrix_add_elem_ui(sparse_M, 2, 3, 1);

    sparse_matrix_add_elem_ui(sparse_M, 3, 3, 6);
    sparse_matrix_add_elem_ui(sparse_M, 3, 4, 1);

    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    printf("gauss rref:\n");
    sparse_matrix_gauss_rref(sparse_M);

    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    sparse_matrix_clear(sparse_M);

    // 0 1 0 0 1 0 0 0 0 0 0 
    // 0 1 6 0 0 0 0 0 0 0 0 
    // 0 0 0 0 0 6 0 0 0 1 0 
    // 0 0 0 0 0 1 0 6 0 0 0 
    // 1 0 0 1 0 0 0 0 0 0 0 
    // 6 0 1 0 0 0 0 0 0 0 0 
    // 0 0 1 0 0 0 1 0 0 0 0 
    // 0 0 0 6 0 0 1 0 0 0 0 
    // 0 0 0 0 6 0 0 0 1 0 0 
    // 0 0 0 0 0 0 6 0 0 0 1 
    // 0 0 0 0 0 0 0 0 1 0 6 

    sparse_matrix_init(sparse_M, 11, 11, field);

    sparse_matrix_add_elem_ui(sparse_M, 0, 1, 1);
    sparse_matrix_add_elem_ui(sparse_M, 0, 4, 1);
    sparse_matrix_add_elem_ui(sparse_M, 1, 1, 1);
    sparse_matrix_add_elem_ui(sparse_M, 1, 2, 6);
    sparse_matrix_add_elem_ui(sparse_M, 2, 5, 6);
    sparse_matrix_add_elem_ui(sparse_M, 2, 9, 1);
    sparse_matrix_add_elem_ui(sparse_M, 3, 5, 1);
    sparse_matrix_add_elem_ui(sparse_M, 4, 0, 1);
    sparse_matrix_add_elem_ui(sparse_M, 4, 3, 1);
    sparse_matrix_add_elem_ui(sparse_M, 5, 2, 1);
    sparse_matrix_add_elem_ui(sparse_M, 6, 2, 1);
    sparse_matrix_add_elem_ui(sparse_M, 3, 7, 6);
    sparse_matrix_add_elem_ui(sparse_M, 5, 0, 6);
    sparse_matrix_add_elem_ui(sparse_M, 6, 6, 1);
    sparse_matrix_add_elem_ui(sparse_M, 7, 3, 6);
    sparse_matrix_add_elem_ui(sparse_M, 7, 6, 1);
    sparse_matrix_add_elem_ui(sparse_M, 8, 4, 6);
    sparse_matrix_add_elem_ui(sparse_M, 8, 8, 1);
    sparse_matrix_add_elem_ui(sparse_M, 9, 6, 6);
    sparse_matrix_add_elem_ui(sparse_M, 9, 10, 1);
    sparse_matrix_add_elem_ui(sparse_M, 10, 8, 1);
    sparse_matrix_add_elem_ui(sparse_M, 10, 10, 6);

    sparse_matrix_print_pretty(sparse_M);

    printf("gauss rref:\n");
    sparse_matrix_gauss_rref(sparse_M);

    sparse_matrix_print_info(sparse_M);
    printf("\n");
    sparse_matrix_print_pretty(sparse_M);
    printf("\n");

    // sparse_matrix_clear(sparse_M);
    // sparse_matrix_init(sparse_M, 19, 27, field);

    // sparse_matrix_add_elem_ui(sparse_M, 0, 2, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 0, 9, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 0, 15, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 1, 2, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 1, 3, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 1, 4, 4);
    // sparse_matrix_add_elem_ui(sparse_M, 1, 6, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 1, 10, 6);
    // sparse_matrix_add_elem_ui(sparse_M, 2, 7, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 2, 17, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 3, 7, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 3, 12, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 3, 14, 4);
    // sparse_matrix_add_elem_ui(sparse_M, 3, 20, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 3, 25, 6);
    // sparse_matrix_add_elem_ui(sparse_M, 4, 5, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 4, 23, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 4, 24, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 5, 5, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 5, 11, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 5, 12, 4);
    // sparse_matrix_add_elem_ui(sparse_M, 5, 14, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 5, 19, 6);
    // sparse_matrix_add_elem_ui(sparse_M, 6, 0, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 6, 5, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 6, 12, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 7, 0, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 7, 7, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 7, 8, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 8, 4, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 8, 13, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 9, 4, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 9, 21, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 9, 22, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 10, 1, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 10, 7, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 10, 14, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 11, 1, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 11, 5, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 12, 3, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 12, 15, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 12, 16, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 13, 4, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 13, 13, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 14, 5, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 14, 23, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 14, 24, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 15, 6, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 15, 15, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 16, 7, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 16, 17, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 17, 8, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 17, 18, 5);
    // sparse_matrix_add_elem_ui(sparse_M, 18, 9, 3);
    // sparse_matrix_add_elem_ui(sparse_M, 18, 13, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 18, 15, 4);
    // sparse_matrix_add_elem_ui(sparse_M, 18, 21, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 18, 26, 6);

    

    // // sparse_matrix_print_pretty(sparse_M);

    // // printf("gauss ref:\n");
    // // r = sparse_matrix_gauss_rref(sparse_M);
    // // printf("rang = %ld\n", r);

    // // sparse_matrix_print_info(sparse_M);
    // // printf("\n");
    // // sparse_matrix_print_pretty(sparse_M);
    // // printf("\n");

    // sparse_matrix_clear(sparse_M);
    // sparse_matrix_init(sparse_M, 2, 5, field);
    
    // sparse_matrix_add_elem_ui(sparse_M, 0, 0, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 0, 2, 2);
    // sparse_matrix_add_elem_ui(sparse_M, 0, 4, 3);

    // sparse_matrix_add_elem_ui(sparse_M, 1, 1, 1);
    // sparse_matrix_add_elem_ui(sparse_M, 1, 3, 1);

    // sparse_matrix_print_pretty(sparse_M);
    // printf("\n");

    // sparse_matrix_add_line_mul_ui(sparse_M, 0, 1, 7);

    // sparse_matrix_print_pretty(sparse_M);
    // printf("\n");

    // sparse_matrix_print_info(sparse_M);
    // printf("\n");

    sparse_matrix_clear(sparse_M);
    fq_nmod_ctx_clear(field);

}