#include "sparse_matrix.h"
#include "types.h"

struct sparse_matrix_pair{
    ulong k;
    fq_nmod_struct* val;
};

typedef struct sparse_matrix_pair sparse_matrix_pair;

void sparse_matrix_pair_init(sparse_matrix_pair* smp, const Field ctx){
    smp->val = flint_calloc(1, sizeof(fq_nmod_struct));
    fq_nmod_init(smp->val, ctx);
}

void sprase_matrix_pair_clear(sparse_matrix_pair* smp, const Field ctx){
    fq_nmod_clear(smp->val, ctx);
    flint_free(smp->val);
}

slong __find_ind(ulong* mas, ulong size, ulong val){
    slong res = -1;

    for(ulong i = 0; i < size; i++)
        if (mas[i] == val){
            res = i;
            break;
        }

    return res;
}

void sparse_matrix_init(sparse_matrix_struct* m, ulong lines, ulong columns, Field ctx){

    m->lines = lines;
    m->columns = columns;  
    m->ctx = ctx; 

    m->l_ind = flint_calloc(lines, sizeof(ulong));
    m->c_ind = flint_calloc(columns, sizeof(ulong));

    ulong i;
    for(i = 0; i < lines; i++)
        m->l_ind[i] = i;

    for(i = 0; i < columns; i++)
        m->c_ind[i] = i;

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
        m->sm[i] = g_array_new(FALSE, FALSE, sizeof(sparse_matrix_pair));

}

void sparse_matrix_clear(sparse_matrix_struct* m){
    sparse_matrix_pair* mas;
    ulong i, j;
    for(i = 0; i < m->size; i++){
        mas = (sparse_matrix_pair*)m->sm[i]->data;
        for(j = 0; j < m->sm[i]->len; j++) sprase_matrix_pair_clear(&mas[j], m->ctx);
        g_array_free(m->sm[i], TRUE);
        
    }

    
    flint_free(m->sm);
    flint_free(m->l_ind);
    flint_free(m->c_ind);
}

// No element availability check
void sparse_matrix_add_elem_fq_nmod(sparse_matrix_struct* m, ulong line, ulong column, fq_nmod_struct* val){
    // Pair p;
    // p.second = val;
    // slong k = __find_ind(m->l_ind, m->lines, line);
    // slong l = __find_ind(m->c_ind, m->columns, column);

    // if (k == -1 || l == -1) return;

    // // printf("column=%ld ,line=%ld , i=%ld, j=%ld\n", column, line, k, l);

    // if (m->main == 0){
    //     p.first = l; 
    //     g_array_append_val(m->sm[k], p);
    // } else {
    //     p.first = k;
    //     g_array_append_val(m->sm[l], p);
    // }

    // m->canonized = 0;
}

// No element availability check
void sparse_matrix_add_elem_ui(sparse_matrix_struct* m, ulong line, ulong column, ulong val){
    sparse_matrix_pair p;
    sparse_matrix_pair_init(&p, m->ctx);
    fq_nmod_set_ui(p.val, val, m->ctx);
    
    slong k = __find_ind(m->l_ind, m->lines, line);
    slong l = __find_ind(m->c_ind, m->columns, column);

    if (k == -1 || l == -1) return;

    // printf("column=%ld ,line=%ld , i=%ld, j=%ld\n", column, line, k, l);

    if (m->main == 0){
        p.k = l; 
        g_array_append_val(m->sm[k], p);
    } else {
        p.k = k;
        g_array_append_val(m->sm[l], p);
    }

    m->canonized = 0;
    // sprase_matrix_pair_clear(&p, m->ctx);
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
    
    for(ulong i = 0; i < m->size; i++){
        el = (Pair*)m->sm[i]->data;
        for(ulong j = 0; j < m->sm[i]->len; j++){
            if (m->main == 0) printf("(%ld, %ld, %ld)\n", i, el->first, el->second);
            else printf("(%ld, %ld, %ld)\n", el->first, i, el->second);
            el++;
        }
    }
    

}


void __sparse_matrix_mas_canonize(sparse_matrix_struct* m, Pair* mas, slong size){
    if (size == 0) return;

    ulong* ind;
    if (m->main == 0) ind = m->c_ind;
    else ind = m->l_ind;

    slong i = 0;
    slong j = size - 1;

    Pair mid = mas[size / 2];
    // printf("i=%ld j=%ld\n", i, j);

    do {
        // printf("%ld %ld\n", mas[i].first, mid.first);
        // printf("%ld %ld\n", ind[mas[i].first], ind[mid.first]);
        while(ind[mas[i].first] < ind[mid.first]) i++;
        while(ind[mas[j].first] > ind[mid.first]) j--;

        // printf("final: %ld %ld\n", ind[mas[i].first], ind[mid.first]);
        // printf("i=%ld j=%ld\n\n", i, j);

        if (i <= j) {
            ulong first = mas[i].first;
            ulong second = mas[i].second;

            mas[i].first = mas[j].first;
            mas[i].second = mas[j].second;

            mas[j].first = first;
            mas[j].second = second;

            i++;
            j--;
        }
    } while (i <= j);

    
    if (j > 0) __sparse_matrix_mas_canonize(m, mas, j + 1);
    if (i < size) __sparse_matrix_mas_canonize(m, &mas[i], size - i);
}

int __is_already_mas_canonized(sparse_matrix_struct* m, Pair* mas, ulong size){
    if (size == 0) return 1;

    ulong* ind;
    if (m->main == 0) ind = m->c_ind;
    else ind = m->l_ind;

    for(ulong i = 0; i < size-1; i++){
        if (ind[mas[i].first] > ind[mas[i+1].first]) return 0;
    }

    return 1;

}

void sparse_matrix_canonize(sparse_matrix_struct* m){
    for(ulong k = 0; k < m->size; k++){
        // printf("line: %ld\n", k);
        // printf("already canonized: %d\n", __is_already_mas_canonized(m, (Pair*)m->sm[k]->data, m->sm[k]->len));
        if (__is_already_mas_canonized(m, (Pair*)m->sm[k]->data, m->sm[k]->len) == 1) continue;
        __sparse_matrix_mas_canonize(m, (Pair*)(m->sm[k]->data), m->sm[k]->len);
    }
    m->canonized = 1;
}

void sparse_matrix_print_pretty(sparse_matrix_struct* m){
    if (m->canonized == 0)
        sparse_matrix_canonize(m);
    
    // printf("LOL: %d\n", m->canonized);

    ulong k;
    Pair* p;

    if (m->main == 0){
        for(ulong i = 0 ; i < m->size; i++){
            k = 0;
            p = (Pair*)m->sm[m->l_ind[i]]->data;
            for(ulong j = 0; j < m->sm[m->l_ind[i]]->len; j++){
                while (k < m->c_ind[p->first]){
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
        ulong* curr_ind = flint_calloc(m->columns, sizeof(ulong));
        ulong i, j;
        ulong k, l;
        Pair p;

        for(i = 0; i < m->columns; i++) curr_ind[i] = 0;

        for(i = 0; i < m->lines; i++){
            k = m->l_ind[i];
            for(j = 0; j < m->columns; j++){
                l = m->c_ind[j];
                
                if (m->sm[l]->len == curr_ind[l]){
                    printf("0 ");
                    continue;
                }
                else{
                    p = g_array_index(m->sm[l], Pair, curr_ind[l]);
                    if (m->l_ind[p.first] > i) {
                        printf("0 ");
                        continue;
                    }
                    printf("%ld ", p.second);
                    ++curr_ind[l];
                }
            }
            printf("\n");
        }
        flint_free(curr_ind);
    }
}

void sparse_matrix_print_info(const sparse_matrix_struct* m){
    printf("canonized=%d\n", m->canonized);
    printf("lines=%ld | columns=%ld | elems=%ld\n", m->lines, m->columns, sparse_matrix_nnz(m));
    printf("lines=[");

    for(ulong i = 0; i < m->lines-1; i++)
        printf("%ld, ", m->l_ind[i]);
    printf("%ld]\n", m->l_ind[m->lines - 1]);

    printf("columns=[");
    for(ulong i = 0; i < m->columns-1; i++)
        printf("%ld, ", m->c_ind[i]);
    printf("%ld]", m->c_ind[m->columns - 1]);
}

void sparse_matrix_swap_columns(sparse_matrix_struct* m, ulong first, ulong second){

    if (first >= m->columns || second >= m->columns) return;

    ulong k = m->c_ind[first];
    m->c_ind[first] = m->c_ind[second];
    m->c_ind[second] = k;

    m->canonized = 0;
}

void sparse_matrix_swap_lines(sparse_matrix_struct* m, ulong first, ulong second){
    if (first >= m->lines || second >= m->lines) return;

    ulong k = m->l_ind[first];
    m->l_ind[first] = m->l_ind[second];
    m->l_ind[second] = k;

    m->canonized = 0;
}

//TODO
void sparse_matrix_rem_item(sparse_matrix_struct* m, ulong i, ulong j){

}

