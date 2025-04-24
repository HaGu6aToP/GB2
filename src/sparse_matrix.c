#include "sparse_matrix.h"
#include "types.h"



void sparse_matrix_init(sparse_matrix_struct* m, ulong lines, ulong columns){

    m->lines = lines;
    m->columns = columns;

    if (lines <= columns){
        m->size = lines;
        m->main = 0;
    } else {
        m->size = columns;
        m->main = 1;
    }
    
    m->canonized = 0;
    m->sm = flint_calloc(m->size, sizeof(GArray*));
    
    for(ulong i = 0; i < m->size; i++)
        m->sm[i] = g_array_new(FALSE, FALSE, sizeof(Pair));

}


void sparse_matrix_clear(sparse_matrix_struct* m){
    for(ulong i = 0; i < m->size; i++){
        g_array_free(m->sm[i], TRUE);
    }
    flint_free(m->sm);
}

void sparse_matrix_add_elem(sparse_matrix_struct* m, ulong line, ulong column, ulong data){
    Pair p;
    p.second = data;

    if (m->main == 0){
        p.first = column;
        g_array_append_val(m->sm[line], p);

    } else {
        p.first = line;
        g_array_append_val(m->sm[column], p);
    }

    m->canonized = 0;
    
}

ulong sparse_matrix_nnz(const sparse_matrix_struct* m){
    ulong res = 0;
    for(ulong k = 0; k < m->size; k++){
        res += m->sm[k]->len;
    }
    return res;
}

void sparse_matrix_print(const sparse_matrix_struct* m){

    printf("lines=%ld | columns=%ld | elems=%ld\n", m->lines, m->columns, sparse_matrix_nnz(m));

    if (m->size == 0){
        printf("{}");
        return;
    }

    Pair* el;
    if (m->main == 0){
        for(ulong i = 0; i < m->size; i++){
            el = (Pair*)m->sm[i]->data;
            for(ulong j = 0; j < m->sm[i]->len; j++){
                printf("(%ld, %ld, %ld)\n", i, el->first, el->second);
                el++;
            }
        }
    }

}

gint cmpCanonize(gconstpointer a, gconstpointer b){
    Pair *A = (Pair*)a, *B = (Pair*)b;
    
    if (A->first < B->first) return -1;
    if (A->first == B->first) return 0;
    else return 1;

    return 0;
}

void sparse_matrix_canonize(sparse_matrix_struct* m){
    for(ulong k = 0; k < m->size; k++)
        g_array_sort(m->sm[k], cmpCanonize);
    m->canonized = 1;
}

// DO: column-main output
void sparse_matrix_print_pretty(sparse_matrix_struct* m){
    if (m->canonized == 0)
        sparse_matrix_canonize(m);
    
    ulong k;
    Pair* p;

    if (m->main == 0){
        for(ulong i = 0 ; i < m->size; i++){
            k = 0;
            p = (Pair*)m->sm[i]->data;
            for(ulong j = 0; j < m->sm[i]->len; j++){
                while (k < p->first){
                    printf("0 ");
                    k++;
                }
                printf("%ld ", p->second);
                p++;
                k++;
            }

            while (k < m->columns){
                printf("0 ");
                k++;
            }

            printf("\n");
        }
    } else {
        sparse_matrix_print(m);
    }
}

void sparse_matrix_rem_item(sparse_matrix_struct* m, ulong i, ulong j){

}

