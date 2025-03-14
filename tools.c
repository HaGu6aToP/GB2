#include "headers.h"

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
    

    for(int i = 0; i < npoli; i++){
        fgets(buff, BUFFER_SIZE, file);
        fq_nmod_mpoly_set_str_pretty(basis[i], buff, variables, ctx);
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