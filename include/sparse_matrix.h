#include "headers.h"
#include "types.h"

struct Entry{
    ulong i;
    ulong j;
    ulong data;
};

struct sparse_matrix_struct{
    GArray** sm;
    ulong lines;
    ulong columns;
    ulong size;
    ulong* l_ind;
    ulong* c_ind;
    int main; // 0 - lines, 1 - columns
    int canonized;
    Field ctx;
};


typedef struct Entry Entry;

typedef struct sparse_matrix_struct sparse_matrix_struct;
typedef sparse_matrix_struct sparse_matrix_t[1];


void sparse_matrix_init(sparse_matrix_struct* m, ulong lines, ulong columns, Field ctx);
void sparse_matrix_clear(sparse_matrix_struct* m);
void sparse_matrix_add_elem_fq_nmod(sparse_matrix_struct* m, ulong line, ulong column, fq_nmod_struct* val);
void sparse_matrix_add_elem_ui(sparse_matrix_struct* m, ulong line, ulong column, ulong val);
void sparse_matrix_print(const sparse_matrix_struct* m);
void sparse_matrix_rem_item(sparse_matrix_struct* m, ulong i, ulong j);
void sparse_matrix_canonize(sparse_matrix_struct* m);
void sparse_matrix_print_pretty(sparse_matrix_struct* m);
void sparse_matrix_swap_columns(sparse_matrix_struct* m, ulong first, ulong second);
void sparse_matrix_swap_lines(sparse_matrix_struct* m, ulong first, ulong second);
void sparse_matrix_add_line_mul_ui(sparse_matrix_struct* m, ulong line, ulong added_line, ulong coeff);
void sparse_matrix_add_line_mul_fq_nmod(sparse_matrix_struct* m, ulong line, ulong added_line, fq_nmod_struct coeff);
void sparse_matrix_mul_line_ui(sparse_matrix_struct* m, ulong i, ulong coeff);
void sparse_matrix_mul_line_fq_nmod(sparse_matrix_struct* m, ulong i, fq_nmod_struct* coeff);
void sparse_matrix_print_info(const sparse_matrix_struct* m);
ulong sparse_matrix_gauss_retucion(sparse_matrix_struct* m);