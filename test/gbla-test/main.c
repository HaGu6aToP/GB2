#include "gbla/gbla_config.h"
#include "gbla/elimination.h"
#include "gbla/matrix.h"
#include "gbla/mapping.h"
#include <stdio.h>

// #include <math.h>

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

void write_jcf_matrix_to_pbm(sm_t *M, const char *fn, int verbose) {
	char buffer[512];
	unsigned char out_byte  = 0;

	ri_t m = M->nrows;
	ci_t n = M->ncols;

	FILE *fh  = fopen(fn, "wb");

	/*  magic PBM header */
#ifdef __LP64__ /*  64bit machine */
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#else /*  32bit machine */
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), fh);

	ri_t i;
	ci_t j, k;
	/*  row width: number of nonzero elements in current row */
	ci_t sz;

	for (i = 0; i < m; ++i) {
		k   = 0;
		sz  = M->rwidth[i];
		for (j = 0; j < n; ++j) {
			if (k < sz && M->pos[i][k] == j) {
				out_byte  |=  (1 << (7 - (j % 8)));
				k++;
			} else {
				out_byte  &=  ~(1 << (7 - (j % 8)));
			}
			if (j % 8 == 7) {
				fwrite(&out_byte, sizeof(unsigned char), 1, fh);
				out_byte  = 0;
			}
		}
		if (j % 8 != 0)
			fwrite(&out_byte, sizeof(unsigned char), 1, fh);

		fflush(fh);
	}
	fclose(fh);
}

sm_t *load_schreyer_matrix(const char *fn, int verbose)
{
  // meta information of matrix
  double density;
  double fs;
  char *fsu;

  // start loading the matrix
  ri_t m;
  ci_t n;
  mod_t     mod;
  ci_t      width;
  int64_t   fl;

  // open in binary mode first to get file size with fseek
  FILE *fh        = fopen(fn,"rb");
  if (fh == NULL) {
    if (verbose > 0)
      printf("File not found!\n");
    return NULL;
  } else {
    fseek(fh, 0L, SEEK_END);
    fl  = ftell(fh);
    fclose(fh);
  }

  // now read data from file
  fh  = fopen(fn,"r");
  // get characteristic
  if (fscanf(fh, "%u", &mod) == 0)
    return NULL;
  // get columns
  if (fscanf(fh, "%u", &m) == 0)
    return NULL;
  // get rows
  if (fscanf(fh, "%u", &n) == 0)
    return NULL;

  // set modulo by hand
  //mod = (mod_t)12451;
  //
  //

  // read entries from file
  sm_t *M   = (sm_t *)malloc(sizeof(sm_t));
  M->rows   = (re_t **)malloc(m*sizeof(re_t *));
  M->pos    = (ci_t **)malloc(m*sizeof(ci_t *));
  M->rwidth = (ci_t *)malloc(m*sizeof(ci_t));

  ri_t i;
  ci_t j;
  ci_t pos;
  re_l_t elt;
  uint64_t nonzeroes =  0;

  for (i = 0; i < m; ++i) {
    // get row width
    if (fscanf(fh, "%u", &width) == 0)
      return NULL;
    M->rwidth[i]  = width;
    // reserve memory in matrix M for rows[i]
    M->rows[i]  = (re_t *)malloc(width * sizeof(re_t));
    M->pos[i]   = (ci_t *)malloc(width * sizeof(ci_t));
    for (j = 0; j < width; ++j) {
      if (fscanf(fh,"%u",&pos) == 0)
        return NULL;
      M->pos[i][j] = pos;
    }
    for (j = 0; j < width; ++j) {
      if (fscanf(fh,"%lu",&elt) == 0)
        return NULL;
      M->rows[i][j] = (re_t)elt;
    }
    nonzeroes   +=  width;
  }
  //
  // density of matrix
  density =   (double) n * (double) m;
  density =   (double) (nonzeroes) / density;
  density *=  100.0;
  // file size of matrix
  fs  = (double) fl / 1024 / 1024;
  fsu = "MB";
  if (fs > 1000) {
    fs  = fs / 1024;
    fsu = "GB";
  }

  // get meta data
  M->nrows    = m;
  M->ncols    = n;
  M->nnz      = nonzeroes;
  M->mod      = mod;
  M->density  = (float)density;
  M->fs       = (float)fs;
  M->fsu      = fsu;


  fclose(fh);
  return M;
}


void main(){
  sm_t* M = load_schreyer_matrix("4.txt", 0);
  sort_schreyer_matrix(M);
  normalize_schreyer_input_rows(M);
  print_sparse_matrix(M);

  reduce_sparse_matrix(M);

  print_sparse_matrix(M);
}


void test(){

  //   // 1, 0, 0, 0, 0, 0, 0, 1, 1, 4, 0, 0
  //   // 1, 6, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0
  //   // 0, 1, 6, 1, 1, 1, 4, 0, 0, 0, 0, 0
  //   // 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 6, 2
    
  //   ri_t i;
  //   // ri_t m = 4;
  //   // ci_t n = 12;

  //   sm_t* M = (sm_t*)malloc(sizeof(sm_t));
  //   // M->mod = 7;
  //   // M->ncols  = n;
  //   // M->nrows  = m;
  //   // M->rows   = (re_t**)malloc(m*sizeof(re_t*));
  //   // M->pos    = (ci_t**)malloc(m*sizeof(ci_t *));
  //   // M->rwidth = (ci_t*)malloc(m*sizeof(ci_t));

  //   // M->rwidth[0] = (ci_t)4;
  //   // M->rwidth[1] = (ci_t)3;
  //   // M->rwidth[2] = (ci_t)6;
  //   // M->rwidth[3] = (ci_t)3;

  //   // for(i = 0; i < m; i++){
  //   //     // M->rows[i] = (re_t*)malloc(M->rwidth[i] * sizeof(re_t));
  //   //     // M->pos[i]  = (ci_t *)malloc(M->rwidth[i] * sizeof(ci_t));
  //   //     M->rows[i] = (re_t*)malloc(n * sizeof(re_t));
  //   //     M->pos[i]  = (ci_t *)malloc(n * sizeof(ci_t));
  //   // }
    
  //   // printf("hello\n");

  //   // M->nnz = 16;

  //   // M->rows[0][0] = 1;
  //   // M->rows[0][1] = 1;
  //   // M->rows[0][2] = 1;
  //   // M->rows[0][3] = 4;

  //   // M->rows[1][0] = 1;
  //   // M->rows[1][1] = 6;
  //   // M->rows[1][2] = 2;

  //   // M->rows[2][0] = 1;
  //   // M->rows[2][1] = 6;
  //   // M->rows[2][2] = 1;
  //   // M->rows[2][3] = 1;
  //   // M->rows[2][4] = 1;
  //   // M->rows[2][5] = 4;

  //   // M->rows[3][0] = 1;
  //   // M->rows[3][1] = 6;
  //   // M->rows[3][2] = 2;

  //   // M->pos[0][0] = 0;
  //   // M->pos[0][1] = 7;
  //   // M->pos[0][2] = 8;
  //   // M->pos[0][3] = 9;

  //   // M->pos[1][0] = 0;
  //   // M->pos[1][1] = 1;
  //   // M->pos[1][2] = 3;

  //   // M->pos[2][0] = 1;
  //   // M->pos[2][1] = 2;
  //   // M->pos[2][2] = 3;
  //   // M->pos[2][3] = 4;
  //   // M->pos[2][4] = 5;
  //   // M->pos[2][5] = 6;

  //   // M->pos[3][0] = 5;
  //   // M->pos[3][1] = 10;
  //   // M->pos[3][2] = 11;

  //   // 1, 0, 1, 0, 0
  //   // 1, 1, 0, 0, 0
  //   // 0, 1, 0, 1, 0
  //   // 0, 0, 0, 6, 1

  //   ri_t m = 4;
  //   ci_t n = 5;

  //   M->ncols  = n;
  //   M->nrows  = m;
  //   M->rows   = (re_t**)malloc(m*sizeof(re_t*));
  //   M->pos    = (ci_t**)malloc(m*sizeof(ci_t*));
  //   M->rwidth = (ci_t*)malloc(m*sizeof(ci_t));

  //   M->rwidth[0] = 2;
  //   M->rwidth[1] = 2;
  //   M->rwidth[2] = 2;
  //   M->rwidth[3] = 2;

  //   for(i = 0; i < m; i++){
  //       M->rows[i] = (re_t*)malloc(M->rwidth[i] * sizeof(re_t));
  //       M->pos[i]  = (ci_t *)malloc(M->rwidth[i] * sizeof(ci_t));
  //   }

  //   M->rows[0][0] = 1;
  //   M->rows[0][1] = 1;
  //   M->rows[1][0] = 1;
  //   M->rows[1][1] = 1;
  //   M->rows[2][0] = 1;
  //   M->rows[2][1] = 1;
  //   M->rows[3][0] = 6;
  //   M->rows[3][1] = 1;

  //   M->pos[0][0] = 0;
  //   M->pos[0][1] = 2;
  //   M->pos[1][0] = 0;
  //   M->pos[1][1] = 1;
  //   M->pos[2][0] = 1;
  //   M->pos[2][1] = 3;
  //   M->pos[3][0] = 3;
  //   M->pos[3][1] = 4;

  //   M->nnz = 8;
  //   M->density = compute_density(M->nnz, M->nrows, M->ncols);
  //   M->mod = 7;
  //   M->fs = 0;
  //   M->fsu = 0;

  //   // write_jcf_matrix_to_pbm(M, "gbla_test.pbm", 0);
  //   // M = load_schreyer_matrix("1.txt", 1);
  //   // M = sort_schreyer_matrix(M);
  //   // normalize_schreyer_input_rows(M);
    
  //   print_sparse_matrix(M);
  //   print_sparse_matrix_info(M);
  //   reduce_sparse_matrix(M);
  //   // write_jcf_matrix_to_pbm(M, "gbla_test.pbm", 0);
  //   print_sparse_matrix_info(M);
  //   print_sparse_matrix(M);



  //   // printf("%ld %ld\n", M->rows[0], M->pos[0]);
    

  //   // for(i = 0; i < m; i++){
  //   //     for(int j = 0; j < n; j++){
  //   //         printf("%d ", M->rows[i][j]);
  //   //     }
  //   //     printf("\n");
  //   // }

  //   // ri_t	ii = 0 ;
	// // for ( ; ii < M->nrows ; ++ii) {
	// // 	if (M->rows[ii] != NULL)
	// // 	free(M->rows[ii]);
	// // 	if (M->pos[ii] != NULL)
	// // 	free(M->pos[ii]);
	// // }

  //   // free(M->rows);
  //   // free(M->pos);
  //   // free(M->rwidth);
  //   // free(M);
  //   // M = NULL;
} 