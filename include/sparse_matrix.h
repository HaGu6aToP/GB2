#include "headers.h"

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
    int main; // 0 - lines, 1 - columns
    int canonized;
};


typedef struct Entry Entry;

typedef struct sparse_matrix_struct sparse_matrix_struct;
typedef sparse_matrix_struct sparse_matrix_t[1];


void sparse_matrix_init(sparse_matrix_struct* m, ulong lines, ulong columns);
void sparse_matrix_clear(sparse_matrix_struct* m);
void sparse_matrix_add_elem(sparse_matrix_struct* M, ulong line, ulong column, ulong data);
void sparse_matrix_print(const sparse_matrix_struct* m);
void sparse_matrix_rem_item(sparse_matrix_struct* m, ulong i, ulong j);
void sparse_matrix_canonize(sparse_matrix_struct* m);
void sparse_matrix_print_pretty(sparse_matrix_struct* m);