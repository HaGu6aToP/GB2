#include "headers.h"

gint cmpPair(gconstpointer a, gconstpointer b){
    Pair *A = (Pair*)a, *B = (Pair*)b;
    printf("\n---------------cmpPair---------------\n");
    printf("A=(%ld, %ld), B=(%ld, %ld)\n", A->first, A->second, B->first, B->second);
    printf("---------------------------------------\n");
    
    if (A->first < B->first) return -1;
    if (A->first > B->first) return 1;

    if (A->second < B->second) return -1;
    if (A->second > B->second) return 1;

    return 0;
}

void log_B(GArray* B){
    printf("|B|=%d, B=\n", (int)B->len);
    for (Pair* i = (Pair*)B->data; i < (Pair*)B->data + B->len; i++){
        printf("(%ld, %ld)", i->first, i->second);
    }
    printf("\n\n");
}

void main(void){
    GArray* B = g_array_new(FALSE, FALSE, sizeof(Pair));

    Pair p1 = {1, 2};
    Pair p2 = {0, 3};
    Pair p3 = {1, 3};
    Pair search = {1, 2};

    g_array_append_val(B, p1);
    g_array_append_val(B, p2);
    g_array_append_val(B, p3);

    log_B(B);


    int index = -1;
    g_array_sort(B, cmpPair);
    int flag = g_array_binary_search(B, &search, cmpPair, &index);

    printf("Search res: %d; index=%d\n", flag, index);
}