#include "suitesparse/cs.h"
#include <stdio.h>


int main(void){
    // 1, 0, 1, 0, 0
    // 1, 1, 0, 0, 0
    // 0, 1, 0, 1, 0
    // 0, 0, 0, 6, 1
    cs_dl *T, *A;
    T = cs_dl_spalloc(0, 0, 1, 1, 1);
    cs_dl_print(T, 0);

    cs_dl_entry(T, 0, 0, 1);
    cs_dl_entry(T, 0, 2, 1);
    cs_dl_entry(T, 1, 0, 1);
    cs_dl_entry(T, 1, 1, 1);
    cs_dl_entry(T, 2, 1, 1);
    cs_dl_entry(T, 2, 3, 1);
    cs_dl_entry(T, 3, 3, 6);
    cs_dl_entry(T, 3, 4, 1);

    cs_dl_print(T, 0);
    A = cs_dl_compress(T);
    cs_dl_print(A, 0);

    cs_dls* S = cs_dl_sqr(4, A, 0);
    double tol = 1e-12;

    cs_dln* res = cs_dl_lu(A, S, tol);

    printf("%ld\n", res);

    // cs_dl_print(res->U, 0);
    
    cs_dl_nfree(res);
    cs_dl_free(T);
    cs_dl_free(A);
}

int test (void)
{
    cs_dl *T, *A, *Eye, *AT, *C, *D ;
    int64_t i, m ;
    int version [3] ;
    cxsparse_version (version) ;
    printf ("CXSparse v%d.%d.%d\n", version [0], version [1], version [2]) ;

    #ifndef TEST_COVERAGE
    if ((version [0] != CS_VER) || (version [1] != CS_SUBVER) ||
        (version [2] != CS_SUBSUB))
    {
        fprintf (stderr, "version in header does not match library\n") ;
        abort ( ) ;
    }
    #endif

    T = cs_dl_load (stdin) ;               /* load triplet matrix T from stdin */
    printf ("T:\n") ; cs_dl_print (T, 0) ; /* print T */
    A = cs_dl_compress (T) ;               /* A = compressed-column form of T */
    printf ("A:\n") ; cs_dl_print (A, 0) ; /* print A */
    cs_dl_spfree (T) ;                     /* clear T */
    AT = cs_dl_transpose (A, 1) ;          /* AT = A' */
    printf ("AT:\n") ; cs_dl_print (AT, 0) ; /* print AT */
    m = A ? A->m : 0 ;                  /* m = # of rows of A */
    T = cs_dl_spalloc (m, m, m, 1, 1) ;    /* create triplet identity matrix */
    for (i = 0 ; i < m ; i++) cs_dl_entry (T, i, i, 1) ;
    Eye = cs_dl_compress (T) ;             /* Eye = speye (m) */
    cs_dl_spfree (T) ;
    C = cs_dl_multiply (A, AT) ;           /* C = A*A' */
    D = cs_dl_add (C, Eye, 1, cs_dl_norm (C)) ;   /* D = C + Eye*norm (C,1) */
    printf ("D:\n") ; cs_dl_print (D, 0) ; /* print D */
    cs_dl_spfree (A) ;                     /* clear A AT C D Eye */
    cs_dl_spfree (AT) ;
    cs_dl_spfree (C) ;
    cs_dl_spfree (D) ;
    cs_dl_spfree (Eye) ;
    return (0) ;
}