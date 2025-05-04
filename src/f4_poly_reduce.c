#include "f4.h" 

#include "gbla/elimination.h"
#include "gbla/matrix.h"
#include "gbla/mapping.h"

void print_sparse_matrix_info(const sm_t* M){
    printf("nrows=%d ncols=%d nnz=%ld density=%f\n", (int)M->nrows, (int)M->ncols, M->nnz, M->density);
    printf("rwidth: ");
    for(int i = 0; i < M->nrows; i++)
        printf("%d ", M->rwidth[i]);
    printf("\n");
}

void print_sparse_matrix(const sm_t* M){
    ri_t i;
    ci_t j;
    int k;
    for(i = 0; i < M->nrows; i++){
        j = 0;
        for(k = 0; k < M->ncols; k++){
            if (j == M->rwidth[i]) break;
            if (M->pos[i][j] == -1){
                printf("0 ");
                ++j;
            } else if (M->pos[i][j] > k){
                printf("0 ");
            } else {
                printf("%d ", M->rows[i][j]);
                ++j;
            }
        }

        while(k < M->ncols){
            printf("0 ");
            ++k;
        }

        printf("\n");
    }

    printf("----\n");

    for(i = 0; i < M->nrows; i++){
        for(j = 0; j < M->rwidth[i]; j++){
            printf("%d ", M->rows[i][j]);
        }
        printf("\n");
    }

    printf("----\n");

    for(i = 0; i < M->nrows; i++){
        for(j = 0; j < M->rwidth[i]; j++){
            printf("%d ", M->pos[i][j]);
        }
        printf("\n");
    }
}

void print_sm_fl_t(const sm_fl_t* M){
    ri_t i;
    ci_t j;
    int k;
    for(i = 0; i < M->nrows; i++){
        j = 0;
        for(k = 0; k < M->ncols; k++){
            if (j == M->sz[i]) break;
            if (M->pos[i][j] == -1){
                printf("0 ");
                ++j;
            } else if (M->pos[i][j] > k){
                printf("0 ");
            } else {
                printf("%d ", M->row[i][j]);
                ++j;
            }
        }

        while(k < M->ncols){
            printf("0 ");
            ++k;
        }

        printf("\n");
    }

    printf("----\n");

    for(i = 0; i < M->nrows; i++){
        for(j = 0; j < M->sz[i]; j++){
            printf("%d ", M->row[i][j]);
        }
        printf("\n");
    }

    printf("----\n");

    for(i = 0; i < M->nrows; i++){
        for(j = 0; j < M->sz[i]; j++){
            printf("%d ", M->pos[i][j]);
        }
        printf("\n");
    }
}

void reduce_sparse_matrix(sm_t* M){
    // construct splicing of matrix M into A, B, C and D
  sb_fl_t *A      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *B     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  sb_fl_t *C      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD
  map_fl_t *map_D = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings for reduced D

  int free_mem = 1, verbose = 0, nthreads = 1;

  splice_fl_matrix_sparse_dense_2(M, A, B, C, D, map, 0, free_mem, verbose, nthreads); 

  // free elements in M, but keep general meta data for later references
	ri_t ii;
	for (ii=0; ii < M->nrows; ++ii) {
		if (M->rows[ii] != NULL)
		free(M->rows[ii]);
		if (M->pos[ii] != NULL)
		free(M->pos[ii]);
	}
  free(M->rows);
  free(M->pos);

  elim_fl_A_sparse_dense_block(&A, B, M->mod, nthreads);
  elim_fl_C_sparse_dense_block(B, &C, D, 1, M->mod, nthreads);

  // copy block D to dense wide (re_l_t) representation
  dm_t *D_red = copy_block_to_dense_matrix(&D, nthreads);
  D_red->mod  = M->mod;

  // eliminate D_red using a structured Gaussian Elimination process on the rows
  ri_t rank_D = 0;
  if (D_red->nrows > 0) rank_D = elim_fl_dense_D(D_red, nthreads);

  reconstruct_matrix_block_no_multiline(M, A, B, D_red, map, nthreads);
}

void F4_poly_reduce(GArray* F_ref, const GArray* F, const GArray* F_monoms, const Field field, const PolynomRing){
    
}