#include "gbla/src/matrix.h"
#include "gbla/src/elimination.h"
#include <stdio.h>
// #include <math.h>

#define GBLA_USE_INT32

// /**
//  * \brief Sparse matrix structure for reading jcf matrices
//  */

//  typedef struct sm_t {
//     ri_t nrows;     /*!<  number of rows */
//     ci_t ncols;     /*!<  number of columns */
//     nnz_t nnz;      /*!<  number of nonzero entries */
//     float density;  /*!<  density used for adjusting memory allocation during
//                           splicing and generation of ABCD blocks */
//     mod_t mod;      /*!<  modulo/field characteristic */
//     float fs;       /*!<  file size of input matrix */
//     char *fsu;      /*!<  file size unit of input matrix, e.g. GB */
//     re_t **rows;    /*!<  address of row: M->rows[i] gives first
//                           address of ith row */
//     ci_t **pos;     /*!<  position of entry in row: M->pos[i] gives first
//                           address of first position of nonzero entry in row i */
//     ci_t *rwidth;   /*!<  width of row: M->rwidth[i] gives number of nonzero
//                           entries in row i */
//     ci_t *buf;      /*!<  stores buffer of memory allocated for the given row*/
//   } sm_t;

// /**
//  * \brief Sparse matrix structure for Faugère-Lachartre decompositions.
//  * For non multiline implementation using small sparse blocks for the A part.
//  */

//  typedef struct sm_fl_t {
//     ri_t nrows;       /*!<  number of rows */
//     ci_t ncols;       /*!<  number of columns */
//     nnz_t nnz;        /*!<  number of nonzero elements */
//     double density;   /*!<  density of this submatrix */
//     re_t **row;       /*!< row entries */
//     ci_t **pos;       /*!< position in row */
//     ci_t *sz;         /*!< size of row */
//     ci_t *buf;        /*!< memory buffer already allocated */
//   } sm_fl_t;


void print_sparse_matrix_info(const sm_t* M){
    printf("nrows=%d ncols=%d nnz=%ld density=%f\n", (int)M->nrows, (int)M->ncols, M->nnz, M->density);
}

void print_sparse_matrix(const sm_t* M){
    ri_t i;
    ci_t j;
    int k;
    for(i = 0; i < M->nrows; i++){
        k = 0;
        j = 0;
        for(k = 0; k < M->ncols; k++){
            if (M->pos[i][j] > k) printf("0 ");
            else{
                printf("%d ", M->rows[i][j]);
                ++j;
            }
        }
        printf("\n");
    }
}

void reduce_sparse_matrix(sm_t* M){
    sb_fl_t *A      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
    dbm_fl_t *B     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
    sb_fl_t *C      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
    dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
    map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t));

    int nthreads = 1;

    splice_fl_matrix_sparse_dense_2(M, A, B, C, D, map, 0, 1, nthreads);

    ri_t ii;
	for (ii=0; ii < M->nrows; ++ii) {
		if (M->rows[ii] != NULL)
		free(M->rows[ii]);
		if (M->pos[ii] != NULL)
		free(M->pos[ii]);
	}
    free(M->rows);
    free(M->pos);

    // // column loops
    // const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
    // const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
    // const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
    // const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
    // // row loops
    // const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
    // const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
    // const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
    // const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

    // printf("rlA=%d, clA=%d\n", rlA, clA);

    // printf("============= AAAAAAAAAAAAAAAA ===================\n");
    // for (int ii=0; ii<rlA; ++ii) {
    //     for (int jj=0; jj<clA; ++jj) {
    //     if (A->blocks[ii][jj].val != NULL) {
    //         printf("%d .. %d\n", ii, jj);
    //         for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
    //         for (int ll=0; ll<A->blocks[ii][jj].sz[kk]; ++ll) {
    //             printf("%d | %d || ", A->blocks[ii][jj].val[kk][ll], A->blocks[ii][jj].pos[kk][ll]);
    //         }
    //         printf("\n");
    //         }
    //     }
    //     }
    // }
    // printf("==================================\n");

    elim_fl_A_sparse_dense_block(&A, B, M->mod, nthreads);
    elim_fl_C_sparse_dense_block(B, &C, D, M->mod, nthreads);

    dm_t *D_red = copy_block_to_dense_matrix(&D, nthreads, 1);
    D_red->mod  = M->mod;

    ri_t rank_D = 0;

    if (D_red->nrows > 0) rank_D = elim_fl_dense_D(D_red, nthreads);

    reconstruct_matrix_block_no_multiline(M, B, D_red, map, nthreads);
}

void main(){

    // 1, 0, 0, 0, 0, 0, 0, 1, 1, 4, 0, 0
    // 1, 6, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0
    // 0, 1, 6, 1, 1, 1, 4, 0, 0, 0, 0, 0
    // 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 6, 2
    
    ri_t i;
    ri_t m = 4;
    ci_t n = 12;

    sm_t* M = (sm_t*)malloc(sizeof(sm_t));
    M->ncols = n;
    M->nrows = m;
    M->rows = (re_t**)malloc(m*sizeof(re_t*));
    M->pos    = (ci_t**)malloc(m*sizeof(ci_t *));
    M->rwidth = (ci_t*)malloc(m*sizeof(ci_t));

    M->rwidth[0] = (ci_t)4;
    M->rwidth[1] = (ci_t)3;
    M->rwidth[2] = (ci_t)6;
    M->rwidth[3] = (ci_t)3;

    for(i = 0; i < 4; i++){
        M->rows[i] = (re_t*)malloc(M->rwidth[i] * sizeof(re_t));
        M->pos[i]  = (ci_t *)malloc(M->rwidth[i] * sizeof(ci_t));
    }
    
    printf("hello\n");

    M->nnz = 15;

    M->rows[0][0] = 1;
    M->rows[0][1] = 1;
    M->rows[0][2] = 1;
    M->rows[0][3] = 4;

    M->rows[1][0] = 1;
    M->rows[1][1] = 6;
    M->rows[1][2] = 2;

    M->rows[2][0] = 1;
    M->rows[2][1] = 6;
    M->rows[2][2] = 1;
    M->rows[2][3] = 1;
    M->rows[2][4] = 1;
    M->rows[2][5] = 4;

    M->rows[3][0] = 1;
    M->rows[3][1] = 6;
    M->rows[3][2] = 2;

    M->pos[0][0] = 0;
    M->pos[0][1] = 7;
    M->pos[0][2] = 8;
    M->pos[0][3] = 9;

    M->pos[1][0] = 0;
    M->pos[1][1] = 1;
    M->pos[1][2] = 3;

    M->pos[2][0] = 1;
    M->pos[2][1] = 2;
    M->pos[2][2] = 3;
    M->pos[2][3] = 4;
    M->pos[2][4] = 5;
    M->pos[2][5] = 6;

    M->pos[3][0] = 5;
    M->pos[3][1] = 10;
    M->pos[3][2] = 11;

    M->density = compute_density(M->nnz, M->nrows, M->ncols);
    M->mod = 7;

    print_sparse_matrix(M);
    print_sparse_matrix_info(M);
    reduce_sparse_matrix(M);
    print_sparse_matrix_info(M);
    print_sparse_matrix(M);


    ri_t	ii = 0 ;
	for ( ; ii < M->nrows ; ++ii) {
		if (M->rows[ii] != NULL)
		free(M->rows[ii]);
		if (M->pos[ii] != NULL)
		free(M->pos[ii]);
	}

    free(M->rows);
    free(M->pos);
    free(M->rwidth);
    free(M);
    M = NULL;
} 