#include "sparse_matrix.h"
#include "types.h"

void sparse_matrix_pair_init(sparse_matrix_pair* smp, const Field ctx){
    smp->val = flint_calloc(1, sizeof(fq_nmod_struct));
    fq_nmod_init(smp->val, ctx);
}

void sparse_matrix_pair_clear(sparse_matrix_pair* smp, const Field ctx){
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

// temporarily only line-main case
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

    
    // remove later
    m->size = lines;
    m->main = 0;
    m->canonized_ind = flint_calloc(lines, sizeof(int));

    for(i = 0; i < lines; i++) m->canonized_ind[i] = 1;

    // if (lines <= columns){
    //     m->size = lines;
    //     m->main = 0;
    // } else {
    //     m->size = columns;
    //     m->main = 1;
    // }

    
    m->canonized = 1;
    m->sm = flint_calloc(m->size, sizeof(GArray*));
    
    for(ulong i = 0; i < m->size; i++)
        m->sm[i] = g_array_new(FALSE, FALSE, sizeof(sparse_matrix_pair));

}

void sparse_matrix_clear(sparse_matrix_struct* m){
    sparse_matrix_pair* mas;
    ulong i, j;
    for(i = 0; i < m->size; i++){
        if (m->sm[i]->len != 0){
            mas = (sparse_matrix_pair*)m->sm[i]->data;
            for(j = 0; j < m->sm[i]->len; j++) sparse_matrix_pair_clear(&mas[j], m->ctx);
        }
        g_array_free(m->sm[i], TRUE);
        
    }
    
    flint_free(m->sm);
    flint_free(m->l_ind);
    flint_free(m->c_ind);
    flint_free(m->canonized_ind);
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


    fmpz_t temp;
    fmpz_init(temp);
    fq_nmod_get_fmpz(temp, val, m->ctx);
    sparse_matrix_add_elem_ui(m, line, column, fmpz_get_ui(temp));
    fmpz_clear(temp);

    
    
}

// No element availability check
// TODO: add column-main canonized_ind
void sparse_matrix_add_elem_ui(sparse_matrix_struct* m, ulong line, ulong column, ulong val){
    sparse_matrix_pair p;
    sparse_matrix_pair_init(&p, m->ctx);
    fq_nmod_set_ui(p.val, val, m->ctx);

    
    
    // printf("LOL %ld\n", p.val->coeffs[0]);

    // slong k = __find_ind(m->l_ind, m->lines, line);
    // slong l = __find_ind(m->c_ind, m->columns, column);
    slong k = m->l_ind[line];
    ulong l = m->c_ind[column];

    // if (k == -1 || l == -1) return;

    // printf("column=%ld ,line=%ld , i=%ld, j=%ld\n", column, line, k, l);

    if (m->main == 0){
        if (m->sm[m->l_ind[line]]->len != 0){
            if (g_array_index(m->sm[m->l_ind[line]], sparse_matrix_pair, m->sm[m->l_ind[line]]->len-1).k > m->c_ind[column]) 
                m->canonized_ind[m->l_ind[line]] = 0;
        } 
    }

    if (m->main == 0){
        p.k = l; 
        g_array_append_val(m->sm[k], p);
    } else {
        p.k = k;
        g_array_append_val(m->sm[l], p);
    }

    // m->canonized = 0;
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

    sparse_matrix_pair* el;
    
    for(ulong i = 0; i < m->size; i++){
        el = (sparse_matrix_pair*)m->sm[i]->data;
        for(ulong j = 0; j < m->sm[i]->len; j++){
            if (m->main == 0){
                // printf("(%ld, %ld, %ld)\n", i, el->first, el->second);
                // printf("(%ld, %ld, ", m->l_ind[i], m->c_ind[el->k]);
                printf("(%ld, %ld, ", i, el->k);
                fq_nmod_print_pretty(el->val, m->ctx);
                printf(")\n");
            } else {
                // printf("(%ld, %ld, %ld)\n", el->first, i, el->second);
                // printf("(%ld, %ld, ", m->l_ind[el->k], m->c_ind[i]);
                printf("(%ld, %ld, ", el->k, i);
                fq_nmod_print_pretty(el->val, m->ctx);
                printf(")\n");
            } 
            el++;
        }
    }
    

}


void __sparse_matrix_mas_canonize(sparse_matrix_struct* m, sparse_matrix_pair* mas, slong size){
    if (size == 0) return;

    ulong* ind;
    ulong ind_size;
    if (m->main == 0){
        ind = m->c_ind;
        ind_size = m->columns;
    }
    else{
        ind = m->l_ind;
        ind_size = m->lines;
    }

    slong i = 0;
    slong j = size - 1;

    sparse_matrix_pair mid = mas[size / 2];
    // printf("i=%ld j=%ld\n", i, j);

    do {
        // printf("%ld %ld\n", mas[i].first, mid.first);
        // printf("%ld %ld\n", ind[mas[i].first], ind[mid.first]);

        // while(ind[mas[i].k] < ind[mid.k]) i++;
        // while(ind[mas[j].k] > ind[mid.k]) j--;
        while(__find_ind(ind, ind_size, mas[i].k) < __find_ind(ind, ind_size, mid.k)) i++;
        while(__find_ind(ind, ind_size, mas[j].k) > __find_ind(ind, ind_size, mid.k)) j--;


        // printf("final: %ld %ld\n", ind[mas[i].first], ind[mid.first]);
        // printf("i=%ld j=%ld\n\n", i, j);

        if (i <= j) {
            ulong first = mas[i].k;
            fq_nmod_struct* second = mas[i].val;

            mas[i].k = mas[j].k;
            mas[i].val = mas[j].val;

            mas[j].k = first;
            mas[j].val = second;

            i++;
            j--;
        }
    } while (i <= j);

    
    if (j > 0) __sparse_matrix_mas_canonize(m, mas, j + 1);
    if (i < size) __sparse_matrix_mas_canonize(m, &mas[i], size - i);
}

int __is_already_mas_canonized(sparse_matrix_struct* m, sparse_matrix_pair* mas, ulong size){
    if (size == 0) return 1;

    ulong* ind;
    ulong ind_size;
    if (m->main == 0){
        ind = m->c_ind;
        ind_size = m->columns;
    }
    else{
        ind = m->l_ind;
        ind_size = m->lines;
    }

    for(ulong i = 0; i < size-1; i++){
        // if (ind[mas[i].k] > ind[mas[i+1].k]) return 0;
        if (__find_ind(ind, ind_size, mas[i].k) > __find_ind(ind, ind_size, mas[i+1].k)) return 0;
    }

    return 1;

}

void sparse_matrix_canonize(sparse_matrix_struct* m){
    for(ulong k = 0; k < m->size; k++){
        // printf("line: %ld\n", k);
        // printf("already canonized: %d\n", __is_already_mas_canonized(m,  (sparse_matrix_pair*)m->sm[k]->data, m->sm[k]->len));
        // if (__is_already_mas_canonized(m, (sparse_matrix_pair*)m->sm[k]->data, m->sm[k]->len) == 1) continue;
        if (m->canonized == 0 || m->canonized_ind[k] == 0){
            __sparse_matrix_mas_canonize(m, (sparse_matrix_pair*)(m->sm[k]->data), m->sm[k]->len);
            m->canonized_ind[k] = 1;
        }
    }

    m->canonized = 1;

}

void sparse_matrix_print_pretty(sparse_matrix_struct* m){
    // sparse_matrix_print_info(m);
    // printf("\n");
    sparse_matrix_canonize(m);
    
    // printf("LOL: %d\n", m->canonized);

    ulong k;
    sparse_matrix_pair* p;

    if (m->main == 0){
        for(ulong i = 0 ; i < m->size; i++){
            k = 0;
            p = (sparse_matrix_pair*)m->sm[m->l_ind[i]]->data;
            for(ulong j = 0; j < m->sm[m->l_ind[i]]->len; j++){
                //m->c_ind[k] < p->k
                while (k < __find_ind(m->c_ind, m->columns, p->k)){
                    printf("0 ");
                    k++;
                }
                // printf("%ld ", p->val);
                fq_nmod_print_pretty(p->val, m->ctx);
                printf(" ");
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
        sparse_matrix_pair p;

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
                    p = g_array_index(m->sm[l], sparse_matrix_pair, curr_ind[l]);
                    // m->l_ind[p.k] > i
                    // __find_ind(m->l_ind, m->lines, p.k) > i
                    if (__find_ind(m->l_ind, m->lines, p.k) > i) {
                        printf("0 ");
                        continue;
                    }
                    // printf("%ld ", p.val);
                    fq_nmod_print_pretty(p.val, m->ctx);
                    printf(" ");
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
    printf("%ld]\n", m->c_ind[m->columns - 1]);

    printf("canonized_ind=[");
    for(ulong i = 0; i < m->size-1; i++)
        printf("%d, ", m->canonized_ind[i]);
    printf("%d]", m->canonized_ind[m->size - 1]);
}


void sparse_matrix_swap_columns(sparse_matrix_struct* m, ulong first, ulong second){

    if (first >= m->columns || second >= m->columns || first < 0 || second < 0) return;

    ulong k = m->c_ind[first];
    m->c_ind[first] = m->c_ind[second];
    m->c_ind[second] = k;

    // if (m->main == 0) m->canonized = 0;
    m->canonized = 0;
    // if (m->main == 0){
    //     for(ulong i = 0; i < m->size; i++) m->canonized_ind[i] = 0;
    // }
}

// TODO: add column-main canonized_ind
void sparse_matrix_swap_lines(sparse_matrix_struct* m, ulong first, ulong second){
    if (first >= m->lines || second >= m->lines || first < 0 || second < 0) return;

    ulong k = m->l_ind[first];
    m->l_ind[first] = m->l_ind[second];
    m->l_ind[second] = k;

    // if (m->main == 1) m->canonized = 0;
    // m->canonized = 0;
}

//TODO
void sparse_matrix_rem_item(sparse_matrix_struct* m, ulong i, ulong j){

}

// line(a) = line(a) + coeff*added_line(b)
// TODO: add column-main canonized_ind, check column-main case
void sparse_matrix_add_line_mul_ui(sparse_matrix_struct* m, ulong line, ulong added_line, ulong coeff){
    if (line >= m->lines || added_line >= m->lines || line < 0 || added_line < 0) return;

    if (coeff % m->ctx->p == 0) return;

    sparse_matrix_canonize(m);

    if (line == added_line) {
        sparse_matrix_mul_line_ui(m, line, 1 + coeff);
        return;
    }

    fq_nmod_t ca, cb, mul;
    fq_nmod_init(ca, m->ctx);
    fq_nmod_init(cb, m->ctx);
    fq_nmod_init(mul, m->ctx);

    if (m->main == 0){
        ulong i, j, len_b, len_a;
        GArray* a = m->sm[m->l_ind[line]];
        len_a = a->len;
        sparse_matrix_pair smp;
        sparse_matrix_pair* b = (sparse_matrix_pair*)m->sm[m->l_ind[added_line]]->data;
        len_b = m->sm[m->l_ind[added_line]]->len;
        
        // printf("len_b=%ld\n", len_b);
        i = 0;
        j = 0;

        while (i < len_b){
            if (j == len_a){
                // printf("__find=%ld\n", __find_ind(m->c_ind, m->columns, b[i].k));
                fq_nmod_set(cb, b[i].val, m->ctx);
                fq_nmod_mul_ui(mul, cb, coeff, m->ctx);
                sparse_matrix_add_elem_fq_nmod(m, line, __find_ind(m->c_ind, m->columns, b[i].k), mul);
                ++i;
            } else {
                smp = g_array_index(a, sparse_matrix_pair, j);
                if (smp.k == b[i].k){
                    // printf("i=%ld\n", i);
                    
                    fq_nmod_set(cb, b[i].val, m->ctx);
                    // fq_nmod_print_pretty(smp.val, m->ctx);

                    fq_nmod_set(ca, smp.val, m->ctx);

                    fq_nmod_mul_ui(mul, cb, coeff, m->ctx);
                    fq_nmod_add(smp.val, ca, mul, m->ctx);

                    if (fq_nmod_is_zero(smp.val, m->ctx) == 1){
                        sparse_matrix_pair_clear(&smp, m->ctx);
                        g_array_remove_index(a, j);
                        --len_a;
                    } else ++j;

                    ++i;
                } else if (__find_ind(m->c_ind, m->columns, b[i].k) < __find_ind(m->c_ind, m->columns, smp.k)){
                    // printf("i=%ld\n", i);
                    sparse_matrix_pair new_pair;
                    sparse_matrix_pair_init(&new_pair, m->ctx);
                    new_pair.k = b[i].k;

                    fq_nmod_set(cb, b[i].val, m->ctx);
                    fq_nmod_mul_ui(mul, cb, coeff, m->ctx);


                    fq_nmod_set(new_pair.val, mul, m->ctx);
                    g_array_insert_val(a, j, new_pair);

                    // printf("LOL\n");

                    ++j;
                    ++len_a;
                    ++i;
                } else ++j;
            }
        }
        // while(i < len_b){
        //     // printf("i=%ld : %ld; j=%ld : %ld\n", i, len_b, j, len_a);
        //     if (j == len_a){
        //         // sparse_matrix_pair new_pair;
        //         // sparse_matrix_pair_init(&new_pair, m->ctx);
        //         // new_pair.k = b[i].k;

        //         // fq_nmod_set(cb, b[i].val, m->ctx);
        //         // fq_nmod_mul_ui(mul, cb, coeff, m->ctx);

        //         // fq_nmod_set(new_pair.val, mul, m->ctx);
        //         // g_array_append_val(a, new_pair);

        //         // m->canonized_ind[m->l_ind[line]] = 0;
        //         sparse_matrix_add_elem_fq_nmod(m, line, __find_ind(m->c_ind, m->columns, b[i].k), b[i].val);

        //         ++i;
        //     } else {
        //         smp = g_array_index(a, sparse_matrix_pair, j);

        //         // printf("smp.k=%ld b[i].k=%ld\n", m->c_ind[smp.k], m->c_ind[b[i].k]);

        //         if (m->c_ind[b[i].k] == m->c_ind[smp.k]){
        //             fq_nmod_set(cb, b[i].val, m->ctx);
        //             fq_nmod_set(ca, smp.val, m->ctx);

        //             fq_nmod_mul_ui(mul, cb, coeff, m->ctx);
        //             fq_nmod_add(smp.val, ca, mul, m->ctx);

        //             if (fq_nmod_is_zero(smp.val, m->ctx) == 1){
        //                 sparse_matrix_pair_clear(&smp, m->ctx);
        //                 g_array_remove_index(a, j);
        //             } else ++j;

        //             ++i;
        //         } else if (m->c_ind[b[i].k] > m->c_ind[smp.k]){
        //             ++j;
        //         } else {
        //             sparse_matrix_pair new_pair;
        //             sparse_matrix_pair_init(&new_pair, m->ctx);
        //             new_pair.k = b[i].k;

        //             fq_nmod_set(cb, b[i].val, m->ctx);
        //             fq_nmod_mul_ui(mul, cb, coeff, m->ctx);

        //             fq_nmod_set(new_pair.val, mul, m->ctx);
        //             printf("LOL %ld %ld : %ld\n", smp.k, j, len_a);
        //             // g_array_append_val(a, new_pair);
        //             g_array_insert_val(a, j, new_pair);

        //             ++j;
        //             ++len_a;

        //             // m->canonized_ind[m->l_ind[line]] = 0;

        //             ++i;
        //         }
        //     }
        // }

    } else {
        ulong i, j, k;
        sparse_matrix_pair* mas;
        fq_nmod_struct* el;
        int finded_line, finded_added_line;

        for(j = 0; j < m->size; j++){
            mas = (sparse_matrix_pair*)m->sm[j]->data;
            finded_line = 0;
            finded_added_line = 0;
            for(i = 0; i < m->sm[j]->len; i++){
                if (mas[i].k == m->l_ind[line]){
                    finded_line = 1;
                    el = mas[i].val;
                    k = i;
                    // printf("finded line, i=%ld\n", i);
                }

                if (mas[i].k == m->l_ind[added_line]){
                    finded_added_line = 1;
                    fq_nmod_set(cb, mas[i].val, m->ctx);
                    // printf("finded added line, i=%ld\n", i);
                }

                if (finded_line && finded_added_line) break;
            }

            // printf("j=%ld, finded_line=%d, finded_added_line=%d\n", j, finded_line, finded_added_line);

            if (finded_line == 1 && finded_added_line == 1){
                fq_nmod_mul_ui(mul, cb, coeff, m->ctx);
                fq_nmod_add(ca, mul, el, m->ctx);
                fq_nmod_set(el, ca, m->ctx);

                if (fq_nmod_is_zero(el, m->ctx) == 1){
                    fq_nmod_clear(el, m->ctx);
                    flint_free(el);
                    g_array_remove_index(m->sm[j], k);
                }

            } else if (finded_line == 0 && finded_added_line == 1){
                sparse_matrix_pair new_pair;
                sparse_matrix_pair_init(&new_pair, m->ctx);
                // printf("k=%ld, true_k=%ld\n", k, m->l_ind[k]);
                new_pair.k = m->l_ind[line];

                fq_nmod_mul_ui(mul, cb, coeff, m->ctx);

                fq_nmod_set(new_pair.val, mul, m->ctx);
                g_array_append_val(m->sm[j], new_pair);
            }

        }
    } 

    fq_nmod_clear(ca, m->ctx);
    fq_nmod_clear(cb, m->ctx);
    fq_nmod_clear(mul, m->ctx);

    // m->canonized = 0;
}

void sparse_matrix_add_line_mul_fq_nmod(sparse_matrix_struct* m, ulong line, ulong added_line, fq_nmod_struct* coeff){
    fmpz_t temp;
    fmpz_init(temp);
    fq_nmod_get_fmpz(temp, coeff, m->ctx);

    // printf("coeff=");
    // fmpz_print(temp);
    // printf("\n");

    sparse_matrix_add_line_mul_ui(m, line, added_line, fmpz_get_ui(temp));

    fmpz_clear(temp);
}

void sparse_matrix_mul_line_ui(sparse_matrix_struct* m, ulong k, ulong coeff){
    fq_nmod_t x;
    fq_nmod_init(x, m->ctx);
    if (m->main == 0){
        sparse_matrix_pair* line = (sparse_matrix_pair*)m->sm[m->l_ind[k]]->data;

        if (coeff % m->ctx->p == 0){
            for(ulong i = 0; i <  m->sm[m->l_ind[k]]->len; i++) sparse_matrix_pair_clear(&line[i], m->ctx);
            g_array_free(m->sm[m->l_ind[k]], TRUE);
            m->sm[m->l_ind[k]] = g_array_new(FALSE, FALSE, sizeof(sparse_matrix_pair));
            // return;
        }

        for(ulong i = 0; i < m->sm[m->l_ind[k]]->len; i++){
            fq_nmod_mul_ui(x, line[i].val, coeff, m->ctx);
            fq_nmod_set(line[i].val, x, m->ctx);
        }

    } else {
        sparse_matrix_pair smp;
        for(ulong j = 0; j < m->size; j++){
            for(ulong i = 0; i < m->sm[j]->len; i++){
                smp = g_array_index(m->sm[j], sparse_matrix_pair, i);
                // m->l_ind[smp.k] == k;
                if (smp.k == m->l_ind[k]){
                    if (coeff % m->ctx->p == 0){
                        sparse_matrix_pair_clear(&smp, m->ctx);
                        g_array_remove_index(m->sm[j], i);
                    } else {
                        fq_nmod_mul_ui(x, smp.val, coeff, m->ctx);
                        fq_nmod_set(smp.val, x, m->ctx);
                    }
                }
            }
        }
    }
    fq_nmod_clear(x, m->ctx);
}

void sparse_matrix_mul_line_fq_nmod(sparse_matrix_struct* m, ulong k, fq_nmod_struct* coeff){
    fmpz_t temp;
    fmpz_init(temp);
    fq_nmod_get_fmpz(temp, coeff, m->ctx);

    sparse_matrix_mul_line_ui(m, k, fmpz_get_ui(temp));

    fmpz_clear(temp);
}

// TODO column-main case
slong sparse_matrix_find_not_null_line(const sparse_matrix_struct* m, ulong begin){
    if (begin < 0) begin = 0;
    if (begin >= m->lines) return -1;
    if (m->main == 0){
        for(ulong i = begin; i < m->size; i++){
            if (m->sm[m->l_ind[i]]->len != 0) return i;
        }
    } else {

    }

    return -1;
}

// TODO column-main case
slong sparse_matrix_find_not_null_el_in_line(sparse_matrix_struct* m, ulong line, ulong begin){
    if (begin < 0) begin = 0;
    if (begin >= m->columns) return -1;
    if (m->main == 0){
        if (m->canonized = 0) sparse_matrix_canonize(m);
        if (m->sm[m->l_ind[line]]->len == 0) return -1;
        return __find_ind(m->c_ind, m->columns, g_array_index(m->sm[m->l_ind[line]], sparse_matrix_pair, 0).k);
    } else {

    }
}


// TODO: column-main case
ulong sparse_matrix_gauss_ref(sparse_matrix_struct* m){
    
    fq_nmod_t x;
    sparse_matrix_pair smp;
    fq_nmod_init(x, m->ctx);
    ulong r;
    sparse_matrix_canonize(m);

    // printf("matrix: \n");
    // sparse_matrix_print_pretty(m);
    // printf("\n");

    if (m->main == 0){
        ulong i = 0, j = 0;
        slong k, l;

        while(1){
            // printf("i=%ld\n", i)
            // if (i > m->lines){
            //     --i;
            //     break;
            // }
            k = sparse_matrix_find_not_null_line(m, i);
            if (k == -1) break;

            if (k != i){
                sparse_matrix_swap_lines(m, k, i);
                // printf("swap lines %ld and %ld\n", i, k);
                // sparse_matrix_print_pretty(m);
                // printf("\n");
            }

            // if (m->sm[m->l_ind[i]]->len == 0){
            //     k = sparse_matrix_find_not_null_line(m, i);
            //     if (k == -1) break;

            //     if (k != i){
            //         sparse_matrix_swap_lines(m, k, i);
            //         printf("swap lines %ld and %ld\n", i, k);
            //         sparse_matrix_print_pretty(m);
            //         printf("\n");
            //     }
            // }
            
            
            l = sparse_matrix_find_not_null_el_in_line(m, i, i);

            if (l != i){
                sparse_matrix_swap_columns(m, l, i);
                sparse_matrix_canonize(m);
                // printf("swap columns %ld and %ld\n", i, l);
                // sparse_matrix_print_pretty(m);
                // printf("\n");
            }

            

            fq_nmod_inv(x, g_array_index(m->sm[m->l_ind[i]] ,sparse_matrix_pair, 0).val ,m->ctx);

            sparse_matrix_mul_line_fq_nmod(m, i, x);

            // printf("mul line %ld by ", i);
            // fq_nmod_print_pretty(x, m->ctx);
            // printf("\n");

            // sparse_matrix_print_pretty(m);
            // printf("\n");

            for(j = i+1; j < m->lines; j++){
                if (m->sm[m->l_ind[j]]->len == 0){
                    // printf("skip line %ld\n", j);
                    continue;
                }

                smp = g_array_index(m->sm[m->l_ind[j]], sparse_matrix_pair, 0);

                if (smp.k != m->c_ind[i]){
                    // printf("skip line %ld\n", j);
                    continue;
                }

                fq_nmod_neg(x, smp.val, m->ctx);

                // printf("x = ");
                // fq_nmod_print_pretty(x, m->ctx);
                // printf("\n");

                // printf("line %ld plus %ld mul ", j, i);
                // printf("x = ");
                // fq_nmod_print_pretty(x, m->ctx);
                // printf("\n");

                // printf("j=%ld, i=%ld\n", j, i);

                sparse_matrix_add_line_mul_fq_nmod(m, j, i, x);
            }

            // sparse_matrix_print_pretty(m);
            // printf("\n");


            sparse_matrix_canonize(m);
            i++;

            // break;
        }

        r = i;

        // while(1){
            
        //     k = sparse_matrix_find_not_null_line(m, i);
        //     // printf("k = %ld\n", k);

        //     if (k == -1) break;

        //     if (k != i) sparse_matrix_swap_lines(m, k, i);

        //     // if (m->canonized == 0) sparse_matrix_canonize(m);

        //     __sparse_matrix_mas_canonize(m, (sparse_matrix_pair*)m->sm[m->l_ind[i]]->data, m->sm[m->l_ind[i]]->len);
        //     l = sparse_matrix_find_not_null_el_in_line(m, i, i);
        //     // printf("l = %ld\n", l);
        //     if (l != i) sparse_matrix_swap_columns(m, l, i);

        //     // printf("%ld %ld\n", m->l_ind[i], m->c_ind[i]);

        //     fq_nmod_inv(x, g_array_index(m->sm[m->l_ind[i]] ,sparse_matrix_pair, 0).val ,m->ctx);
        //     // printf("x = ");
        //     // fq_nmod_print_pretty(x, m->ctx);
        //     // printf("\n");

        //     sparse_matrix_mul_line_fq_nmod(m, i, x);


        //     sparse_matrix_canonize(m);
        //     // sparse_matrix_print_pretty(m);
        //     // printf("\n");

        //     for(j = i+1; j < m->lines; j++){
        //         if (m->sm[m->l_ind[j]]->len == 0){
        //             // printf("skip line %ld\n", j);
        //             continue;
        //         }
        //         smp = g_array_index(m->sm[m->l_ind[j]], sparse_matrix_pair, 0);
        //         if (smp.k != m->c_ind[i]){
        //             // printf("skip line %ld\n", j);
        //             continue;
        //         }

        //         fq_nmod_neg(x, smp.val, m->ctx);

        //         // printf("x = ");
        //         // fq_nmod_print_pretty(x, m->ctx);
        //         // printf("\n");

        //         // printf("j=%ld, i=%ld\n", j, i);

        //         sparse_matrix_add_line_mul_fq_nmod(m, j, i, x);
        //     }

        //     sparse_matrix_canonize(m);
        //     // sparse_matrix_print_pretty(m);
        //     // printf("\n");

        //     // break;
        //     i++;
        // }
        
        
        // r = i;
    } else {

    }

    
    fq_nmod_clear(x, m->ctx);
    return r;
}

// TODO columns-main case
void sparse_matrix_print_in_file(sparse_matrix_struct* m, char* filename){
    FILE* file;

    if ((file = fopen(filename, "w")) == NULL){
        printf("Failed to open file\n");
        return;
    }

    sparse_matrix_canonize(m);

    ulong k;
    sparse_matrix_pair* p;
    if (m->main == 0){
        for(ulong i = 0 ; i < m->size; i++){
            k = 0;
            p = (sparse_matrix_pair*)m->sm[m->l_ind[i]]->data;
            for(ulong j = 0; j < m->sm[m->l_ind[i]]->len; j++){
                //m->c_ind[k] < p->k
                while (k < __find_ind(m->c_ind, m->columns, p->k)){
                    fprintf(file, "0 ");
                    k++;
                }
                // printf("%ld ", p->val);
                fq_nmod_fprint_pretty(file, p->val, m->ctx);
                fprintf(file, " ");
                p++;
                k++;
            }

            while (k < m->columns){
                fprintf(file, "0 ");
                k++;
            }

            fprintf(file, "\n");
        }
    } else {

    }

    fclose(file);
}

// TODO: columns-main case
ulong sparse_matrix_gauss_rref(sparse_matrix_struct* m){
    ulong r;
    if (m->main == 0){ 
        r = sparse_matrix_gauss_ref(m);
        // sparse_matrix_print_in_file(m, "res.txt");
        // printf("r=%ld\n", r);
        ulong k;
        ulong i;
        fq_nmod_t x;
        fq_nmod_init(x, m->ctx);
        GArray* line;
        sparse_matrix_pair smp;
        for(i = 0; i < r-1; i++){
            // printf("i=%ld\n", i);
            line = m->sm[m->l_ind[i]];

            // sparse_matrix_pair* psmp = (sparse_matrix_pair*)line->data;
            // for(ulong j = 0; j < line->len; j++){
            //     fq_nmod_print_pretty(psmp[j].val, m->ctx);
            //     printf(" ");
            // }

            // printf("\n");

            while(1) {
                if (line->len == 1) break;
                smp = g_array_index(line, sparse_matrix_pair, 1);
                k = __find_ind(m->c_ind, m->columns, smp.k);
                if (k >= r) break;
                fq_nmod_neg(x, smp.val, m->ctx);
                sparse_matrix_add_line_mul_fq_nmod(m, i, k, x);
            }
        }

        fq_nmod_clear(x, m->ctx);
    }

    return r;
}

// TODO: column-main case
ulong sparse_matrix_els_in_line(const sparse_matrix_struct* m, ulong line){
    // if (line < 0 || line >= m->lines) return 0;
    if (m->main == 0) return m->sm[m->l_ind[line]]->len;
}

// TODO: column-main case
void sparse_matrix_true_entry_fq_nmod(sparse_matrix_pair* smp, const sparse_matrix_struct* m, ulong i, ulong position){
    if (m->main == 0){
        // sparse_matrix_pair* mas = (sparse_matrix_pair*)m->sm[m->l_ind[i]]->data;
        // for(ulong j = 0; j < m->sm[m->l_ind[i]]->len; i++){
        //     if (mas[j].k == m->c_ind[j]) return mas[j].val;
        // }
        // return &g_array_index(m->sm[m->l_ind[i]], sparse_matrix_pair, position);
        sparse_matrix_pair temp = g_array_index(m->sm[m->l_ind[i]], sparse_matrix_pair, position);
        smp->k = __find_ind(m->c_ind, m->columns, temp.k);
        smp->val = temp.val;
    }
}