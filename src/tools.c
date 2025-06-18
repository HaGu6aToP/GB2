#include "tools.h"

ulong max(ulong a, ulong b){
    if (a > b) return a;
    else return b;
}

gint cmpPair(gconstpointer a, gconstpointer b){
    Pair *A = (Pair*)a, *B = (Pair*)b;
    // printf("\n---------------cmpPair---------------\n");
    // printf("A=(%ld, %ld), B=(%ld, %ld)\n", A->a, A->b, B->a, B->b);
    // printf("---------------------------------------\n");
    
    if (A->first < B->first) return -1;
    if (A->first > B->first) return 1;

    if (A->second < B->second) return -1;
    if (A->second > B->second) return 1;

    return 0;
}



ulong max_poly_in_lst(const GArray* g, PolynomRing ctx){
    Polynom *max, *f;
    ulong res;

    max = (Polynom*)g->data;
    res = 0;
    f = (Polynom*)g->data;
    f++;
    for(ulong i = 1; i < g->len; i++){
        if (fq_nmod_mpoly_cmp(*f, *max, ctx) == 1){
            res = i;
            max = f;
        }
        f++;
    }
    return res;
}

void *str_key_destroyer(gpointer data){
    free(data);
    return;
}

Polynom max_poly_in_GHashtable(const GHashTable* hash_table, PolynomRing ctx){
    Polynom max, f;
    GHashTableIter i;

    g_hash_table_iter_init(&i, hash_table);
    g_hash_table_iter_next(&i, NULL, &max);


    while(g_hash_table_iter_next(&i, NULL, &f)){
        if (fq_nmod_mpoly_cmp(f, max, ctx) == 1){
            max = f;
        }
    }

    return max;
}


void get_variables(const char** variables, ulong nvars, const char* str){
    ulong var_len = 0;
    char buff[BUFFER_SIZE];
    char* var;
    int counter = 0;
    if (variables == NULL){
        printf("get_variable memory error\n");
        return;
    }

    for(int i = 0; i < BUFFER_SIZE; i++){
        
        if (str[i] == ' ' || str[i] == '\n'){
            buff[var_len] = '\0';
            var = calloc(var_len+1, sizeof(char));
            strcpy(var, buff);
            variables[counter] = var;
            counter++;
            var_len = 0;
        } else {
            buff[var_len] = str[i];
            var_len++;
        }

        if (str[i] == '\n'){
            return;
        }
    }
}

void free_variables(const char** variables, ulong nvars){
    for (int i = 0; i < nvars; i++)
        free((void*)variables[i]);
    free((void*)variables);
}

void read_polinomials(Basis basis, ulong npoli,  const char** variables, PolynomRing ctx, FILE* file){
    char buff[BUFFER_SIZE];
    
    // for(int i = 0; i < 3; i++)
    //     printf("%s-\n", variables[i]);
    for(int i = 0; i < npoli; i++){
        fgets(buff, BUFFER_SIZE, file);
        
        if (buff[strlen(buff) - 1] == '\n') buff[strlen(buff) - 1] = '\0';

        // printf("poli = %s", buff);
        // fq_nmod_mpoly_print_pretty(basis[i], variables, ctx);
        // printf(" %d ", i);


        fq_nmod_mpoly_set_str_pretty(basis[i], buff, variables, ctx);

        // fq_nmod_mpoly_print_pretty(basis[i], variables, ctx);
        // printf("\n");
        // fq_nmod_mpoly_one(basis[i], ctx);
    }
}

void log_B(GArray* B){
    printf("|B|=%d, B=\n", (int)B->len);
    for (Pair* i = (Pair*)B->data; i < (Pair*)B->data + B->len; i++){
        printf("(%ld, %ld)", i->first, i->second);
    }
    printf("\n\n");
}

void log_G(GArray* G, PolynomRing ctx){
    printf("\nG=\n");
    for (Polynom* p = (Polynom*)G->data; p < (Polynom*)G->data + G->len; p++){
        fq_nmod_mpoly_print_pretty(*p, NULL, ctx);
        printf("\n");
    }
    printf("\n");
}

int parseInt(char* chars)
{
    int sum = 0;
    int len = strlen(chars);
    for (int x = 0; x < len; x++)
    {
        int n = chars[len - (x + 1)] - '0';
        sum = sum + powInt(n, x);
    }
    return sum;
}

int powInt(int x, int y)
{
    for (int i = 0; i < y; i++)
    {
        x *= 10;
    }
    return x;
}

void print_poly(const char* header, Polynom p, const char** vars, PolynomRing ctx){
    printf("%s\n", header);
    fq_nmod_mpoly_print_pretty(p, vars, ctx);
    printf("\n");
}

void init_SPair(SPair* pspair, PolynomRing ctx){
    pspair->poly = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    fq_nmod_mpoly_init(pspair->poly, ctx);
}

void set_SPair(SPair* pspair, Polynom p, ulong first, ulong second, PolynomRing ctx){
    fq_nmod_mpoly_set(pspair->poly, p, ctx);
    pspair->first = first;
    pspair->second = second;
}

void free_SPair(SPair* pspair, PolynomRing ctx){
    fq_nmod_mpoly_clear(pspair->poly, ctx);
    flint_free(pspair->poly);
}

void copy_SPair(SPair* pspair, const SPair* resourse){
    pspair->poly = resourse->poly;
    pspair->first = resourse->first;
    pspair->second = resourse->second;
}

void print_ulong_garray(GArray* g){
    for(int i = 0; i < g->len - 1; i++)
        printf("%ld ", g_array_index(g, ulong, i));
    printf("%ld\n", g_array_index(g, ulong, g->len-1));
}

void quick_sort(ulong *s_arr, int first, int last)
{
    if (first < last)
    {
        int left = first, right = last, middle = s_arr[(left + right) / 2];
        do
        {
            while (s_arr[left] < middle) left++;
            while (s_arr[right] > middle) right--;
            if (left <= right)
            {
                int tmp = s_arr[left];
                s_arr[left] = s_arr[right];
                s_arr[right] = tmp;
                left++;
                right--;
            }
        } while (left <= right);
        quick_sort(s_arr, first, right);
        quick_sort(s_arr, left, last);
    }
}

void poly_quick_sort(GArray* g, int first, int last, int rev, const PolynomRing ctx){
    if (first < last)
    {
        int p;
        if (rev == 0)
            p = -1;
        else
            p = 1;
        int left = first, right = last;
        Polynom left_p, right_p;
        Polynom middle = g_array_index(g, Polynom, (left + right) / 2);
        do
        {
            while(fq_nmod_mpoly_cmp(g_array_index(g, Polynom, left), middle, ctx) == p)
                left++;
            while(fq_nmod_mpoly_cmp(g_array_index(g, Polynom, right), middle, ctx) == -p)
                right--;
            // while (s_arr[left] < middle) left++;
            // while (s_arr[right] > middle) right--;
            
            if (left <= right)
            {
                Polynom tmp = g_array_index(g, Polynom, left);
                Polynom* l = ((Polynom*)g->data) + left;
                *l = g_array_index(g, Polynom, right);
                l = ((Polynom*)g->data) + right;
                *l = tmp;
                // int tmp = s_arr[left];
                // s_arr[left] = s_arr[right];
                // s_arr[right] = tmp;
                left++;
                right--;
            }
        } while (left <= right);
        poly_quick_sort(g, first, right, rev, ctx);
        poly_quick_sort(g, left, last, rev, ctx);
    }
}

ulong sum(ulong* arr, ulong len){
    ulong res = 0;
    for(ulong i = 0; i < len; i++){
        res += arr[i];
    }
    return res;
}

slong poly_binary_search(const GArray* g, const Polynom p, const PolynomRing ctx){
    slong l = 0;
    slong r = g->len - 1;
    slong m;
    while( r >= l){
        m = (l + r)/2;
        if (fq_nmod_mpoly_equal(g_array_index(g, Polynom, m), p, ctx) == 1) return m;

        if (fq_nmod_mpoly_cmp(p, g_array_index(g, Polynom, m), ctx) == -1) r = m - 1;
        else l = m + 1;
    }
    return -1;
}

int is_poly_in_lst(const GArray* g, const Polynom p, const PolynomRing ctx){
    Polynom* hp;
    ulong i;
    hp = (Polynom*)g->data;
    for(i = 0; i < g->len; i++){
        if (fq_nmod_mpoly_equal(p, *hp, ctx) == 1)
            return 1;
        hp++;
    }
    return 0;
}

#define mmix(h,k) { k *= m; k ^= k >> r; k *= m; h *= m; h ^= k; }
unsigned int MurmurHash2A ( const void * key, int len, unsigned int seed )
{
	const unsigned int m = 0x5bd1e995;
	const int r = 24;
	unsigned int l = len;

	const unsigned char * data = (const unsigned char *)key;

	unsigned int h = seed;
	unsigned int k;

	while(len >= 4)
	{
		k = *(unsigned int*)data;

		mmix(h,k);

		data += 4;
		len -= 4;
	}

	unsigned int t = 0;

	switch(len)
	{
	case 3: t ^= data[2] << 16;
	case 2: t ^= data[1] << 8;
	case 1: t ^= data[0];
	};

	mmix(h,t);
	mmix(h,l);

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	return h;
}

void monom_lst_from_poly_lst(GArray* res, const GArray* g, const PolynomRing ctx){
    Polynom* hp;
    Polynom h;
    Polynom f;
    // fq_nmod_mpoly_t m;
    // fq_nmod_mpoly_init(m, ctx);

    hp = (Polynom*)g->data;
    char* key;
    
    // printf("hash\n");

   GHashTable* gh = g_hash_table_new_full(g_str_hash, g_direct_equal, str_key_destroyer, NULL);
   for(ulong i = 0; i < g->len; i++){
        for(ulong j = 0; j < fq_nmod_mpoly_length(*hp, ctx); j++){
            h = flint_calloc(1, sizeof(fq_nmod_mpoly_t));
            fq_nmod_mpoly_init(h, ctx);
            fq_nmod_mpoly_get_term_monomial(h, *hp, j, ctx);
            key = fq_nmod_mpoly_get_str_pretty(h, NULL, ctx);
            if (g_hash_table_lookup(gh, key) == NULL){
                g_hash_table_insert(gh, key, h);
                // printf("%s\n", key);
            }
            // free(key);
        }
        hp++;
   }
//    printf("ok\n");

   GPtrArray* vals = g_hash_table_get_values_as_ptr_array(gh);
   g_array_insert_vals(res, 0, vals->pdata, g_hash_table_size(gh));

   g_ptr_array_free(vals, TRUE);
   g_hash_table_destroy(gh);

//    for(int i = 0; i < g_hash_table_size(gh); i++){
//         fq_nmod_mpoly_print_pretty(((Basis)(res->data))[i], NULL, ctx);
//         printf("\n");
//    }

   

//    while(1){}

    // printf("-------------------monom_lst_from_poly_lst-------------------------\n");
    // print_poly_lst(g, ctx);
    // printf("\n");

    // for(int i = 0; i < g->len; i++){
    //     for(int j = 0; j < fq_nmod_mpoly_length(*hp, ctx); j++){
    //         fq_nmod_mpoly_get_term_monomial(m, *hp, j, ctx);
    //         if (is_poly_in_lst(res, m, ctx) == 0){
    //             f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
    //             fq_nmod_mpoly_init(f, ctx);
    //             fq_nmod_mpoly_set(f, m, ctx);
    //             g_array_append_val(res, f);
    //         } else {
    //             printf("already in\n");
    //         }
    //     }
    //     hp++;
    // }

    // fq_nmod_mpoly_clear(m, ctx);
}

void head_monom_lst_from_poly_lst(GArray* res, const GArray* g, const PolynomRing ctx){
    Polynom* hp;
    fq_nmod_mpoly_t m;
    Polynom f;

    fq_nmod_mpoly_init(m, ctx);
    hp = (Polynom*)g->data;

    for(int i = 0; i < g->len; i++){
        fq_nmod_mpoly_get_term_monomial(m, *hp, 0, ctx);
        if (is_poly_in_lst(res, m, ctx) == 0){
            f = flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
            fq_nmod_mpoly_init(f, ctx);
            fq_nmod_mpoly_set(f, m, ctx);
            g_array_append_val(res, f);
        }
        hp++;
    }

    fq_nmod_mpoly_clear(m, ctx);
}

void free_poly_lst(GArray* g, PolynomRing ctx){
    // Polynom f;
    // printf("%d\n", g->len);
    // ulong len = g->len;
    // for(ulong i = 0; i < len; i++){
    //     f = g_array_index(g, Polynom, g->len-1);
    //     // printf("%ld\n", f);
    //     fq_nmod_mpoly_clear(f, ctx);
    //     flint_free(f);
    //     g_array_remove_index(g, g->len-1);
    // }
    // g_array_free(g, TRUE);
    Polynom* pf = (Polynom*)g->data;
    for(ulong i = 0; i < g->len; i++) {
        fq_nmod_mpoly_clear(pf[i], ctx);
        flint_free(pf[i]);
    }
    g_array_free(g, TRUE);
}

void print_poly_lst(const GArray* lst, const PolynomRing ctx){
    if (lst->len == 0){
        printf("{}");
        return;
    }

    printf("{ ");
    for(int i = 0; i < lst->len-1; i++){
        fq_nmod_mpoly_print_pretty(g_array_index(lst, Polynom, i), NULL, ctx);
        printf(",\n");
    }

    fq_nmod_mpoly_print_pretty(g_array_index(lst, Polynom, lst->len - 1), NULL, ctx);
    printf(" }");
}

void print_hash_table(const GHashTable* hash_table, const PolynomRing ctx){
    if (g_hash_table_size(hash_table) == 0){
        printf("{}");
        return;
    }

    GPtrArray* vals = g_hash_table_get_values_as_ptr_array(hash_table);

    printf("{ ");
    for(int i = 0; i < vals->len-1; i++){
        fq_nmod_mpoly_print_pretty((Polynom)vals->pdata[i], NULL, ctx);
        printf(",\n");
    }

    fq_nmod_mpoly_print_pretty((Polynom)vals->pdata[vals->len - 1], NULL, ctx);
    printf(" }");

    g_ptr_array_free(vals, FALSE);
}

void* __calloc_poly_lst(){
    return g_array_new(FALSE, FALSE, sizeof(Polynom));
}
void* __calloc_poly(){
    return flint_calloc(1, sizeof(fq_nmod_mpoly_struct));
}
