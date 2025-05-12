#include "f4.h" 

#include "gbla/elimination.h"
#include "gbla/matrix.h"
#include "gbla/mapping.h"

#define __DEBUG_F4_POLY_REDUCE 0

void print_sparse_matrix_info(const sm_t* M){
    printf("nrows=%d ncols=%d nnz=%ld density=%f\n", (int)M->nrows, (int)M->ncols, M->nnz, M->density);
    // printf("rwidth: ");
    // for(int i = 0; i < M->nrows; i++)
    //     printf("%d ", M->rwidth[i]);
    printf("\n");
}

void print_sparse_matrix_to_file(char* filename, const sm_t* M){
    FILE* f = fopen(filename, "w");
    ri_t i;
    ci_t j;
    int k;
    for(i = 0; i < M->nrows; i++){
        j = 0;
        for(k = 0; k < M->ncols; k++){
            if (j == M->rwidth[i]) break;
            if (M->pos[i][j] == -1){
                fprintf(f, "0 ");
                ++j;
            } else if (M->pos[i][j] > k){
                fprintf(f, "0 ");
            } else {
                fprintf(f, "%d ", (int)M->rows[i][j]);
                ++j;
            }
        }
  
        while(k < M->ncols){
            fprintf(f, "0 ");
            ++k;
        }
  
        fprintf(f, "\n");
    }
    fclose(f);
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

void print_map(const map_fl_t* map, const sm_t* M){
    printf("\n=======map-info=========\n");
    printf("npiv=%d\n", map->npiv);
    printf("pc = ");
    for(int i = 0; i < M->ncols; i++) printf("%d ", map->pc[i]);
    printf("\nnpc= ");
    for(int i = 0; i < M->ncols; i++) printf("%d ", map->npc[i]);
    printf("\npc_rev= ");
    for(int i = 0; i < M->ncols; i++) printf("%d ", map->pc_rev[i]);
    printf("\nnpc_rev= ");
    for(int i = 0; i < M->ncols; i++) printf("%d ", map->npc_rev[i]);
    printf("\npri= ");
    for(int i = 0; i < M->nrows; i++) printf("%d ", map->pri[i]);
    printf("\nnpri= ");
    for(int i = 0; i < M->nrows; i++) printf("%d ", map->npri[i]);
    printf("\n=======map-info=========\n");
  }

ulong* reduce_sparse_matrix(sm_t* M){
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

  ulong* res = calloc(M->ncols, sizeof(ulong));
  ulong i;
  for(i = 0; i < map->npiv; i++) res[i] = map->pc_rev[i];
  for(i; i < M->ncols; i++) res[i] = map->npc_rev[i - map->npiv];

//   #if __DEBUG_F4_POLY_REDUCE
//     printf("map:\n");
//     print_map(map, M);
//     printf("\n");
//   #endif
  

  reconstruct_matrix_block_no_multiline(M, A, B, D_red, map, nthreads);
  return res;
}

void F4_poly_reduce(GArray* F_ref, const GArray* F, const GArray* F_monoms, const Field field, const PolynomRing ctx){
    sm_t* M;
    ulong i, j, k, l;
    fq_nmod_mpoly_t m, sum;
    Basis b, monoms;
    fq_nmod_t coeff;
    Polynom new_poly;
    fmpz_t f;
//-------------------------------------------------------
    M = (sm_t*)malloc(sizeof(sm_t));
    M->mod = fmpz_get_ui(&field->p);
    M->ncols = F_monoms->len;
    M->nrows = F->len;
    M->rows = (re_t**)malloc(M->nrows*sizeof(re_t*));
    M->pos = (ci_t**)malloc(M->nrows*sizeof(ci_t*));
    M->rwidth = (ci_t*)malloc(M->nrows*sizeof(ci_t));

    fq_nmod_mpoly_init(m, ctx);
    fq_nmod_mpoly_init(sum, ctx);
    fq_nmod_init(coeff, field);
    fmpz_init(f);

    b = (Basis)F->data;
    monoms = (Basis)F_monoms->data;
//-------------------------------------------------------
    for(i = 0; i < F->len; i++){
        l = fq_nmod_mpoly_length(b[i], ctx);

        M->rwidth[i] = l;
        M->rows[i] = (re_t*)malloc(l*sizeof(re_t));
        M->pos[i] = (ci_t*)malloc(l*sizeof(ci_t));

        for(j = 0; j < l; j++){
            fq_nmod_mpoly_get_term_monomial(m, b[i], j, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, b[i], j, ctx);
            
            for(k = 0; k < F_monoms->len; k++){
                if (fq_nmod_mpoly_equal(m, monoms[k], ctx)){
                    fq_nmod_get_fmpz(f, coeff, field);
                    M->rows[i][j] = fmpz_get_ui(f);
                    M->pos[i][j] = k;
                    break;
                }
            }
            
        }
    }

    sort_schreyer_matrix(M);
    normalize_schreyer_input_rows(M);

    #if __DEBUG_F4_POLY_REDUCE
        ulong nnz = 0;
        for(ulong i = 0; i < M->nrows; i++)
            nnz += M->rwidth[i];

        M->nnz = nnz;

        printf("M:\n");
        M->density = compute_density(M->nnz, M->nrows, M->ncols);
        print_sparse_matrix_info(M);
        print_sparse_matrix(M);

        if (M->ncols > 500 || M->nrows > 500) print_sparse_matrix_to_file("TT.txt", M);
    #endif

    ulong* p = reduce_sparse_matrix(M);
    
    #if __DEBUG_F4_POLY_REDUCE
    //     printf("columns order:\n");
    //     for(ulong i = 0; i < F_monoms->len; i++) printf("%ld ", p[i]);
        printf("\n reduced M:\n");
        print_sparse_matrix(M);
    #endif


    for(i = 0; i < M->nrows; i++){
        new_poly = flint_calloc(1, sizeof(fq_nmod_mpoly_t));
        fq_nmod_mpoly_init(new_poly, ctx);
        fq_nmod_mpoly_zero(new_poly, ctx);

        for(j = 0; j < M->rwidth[i]; j++){
            // fq_nmod_mpoly_set(m, g_array_index(F_monoms, Polynom, p[M->pos[i][j]]), ctx);
            fq_nmod_set_ui(coeff, M->rows[i][j], field);
            fq_nmod_mpoly_scalar_mul_fq_nmod(m, g_array_index(F_monoms, Polynom, p[M->pos[i][j]]), coeff, ctx);
            fq_nmod_mpoly_add(sum, new_poly, m, ctx);
            fq_nmod_mpoly_set(new_poly, sum, ctx);
        }

        g_array_append_val(F_ref, new_poly);
    }
//-------------------------------------------------------
    fq_nmod_mpoly_clear(m, ctx);
    fq_nmod_mpoly_clear(sum, ctx);
    fq_nmod_clear(coeff, field);
    fmpz_clear(f);

    free(M->rows);
    free(M->pos);
    free(M->rwidth);
    free(M);
    free(p);
}
