#include "mat_utils.h"

/*-------------------------------------------------------------*/
mat_t* init_mat(int nrows, int ncols)
/*-------------------------------------------------------------*/
{
  mat_t *mat = (mat_t*) malloc(sizeof(mat_t));
  mat->nr = nrows;
  mat->nc = ncols;
  mat->m =
    (double**) malloc(sizeof(double*) * mat->nr);
  int i;
  for(i = 0; i < mat->nr; ++i){
    mat->m[i] =
      (double*) malloc(sizeof(double) * mat->nc);
  }
  return mat;
}

/*-------------------------------------------------------------*/
mat_int_t* init_mat_int(int nrows, int ncols)
/*-------------------------------------------------------------*/
{
  mat_int_t *mat = (mat_int_t*) malloc(sizeof(mat_t));
  mat->nr = nrows;
  mat->nc = ncols;
  mat->m =
    (int**) malloc(sizeof(int*) * mat->nr);
  int i;
  for(i = 0; i < mat->nr; ++i){
    mat->m[i] =
      (int*) malloc(sizeof(int) * mat->nc);
  }
  return mat;
}

/*-------------------------------------------------------------*/
void free_mat(mat_t* mat)
/*-------------------------------------------------------------*/
{
  int i;
  for(i = 0; i < mat->nr; ++i){
    free(mat->m[i]);
  }
  free(mat->m);
  free(mat);
  return;
}

/*-------------------------------------------------------------*/
void free_mat_int(mat_int_t* mat)
/*-------------------------------------------------------------*/
{
  int i;
  for(i = 0; i < mat->nr; ++i){
    free(mat->m[i]);
  }
  free(mat->m);
  free(mat);
  return;
}

/*-------------------------------------------------------------*/
void free_mat_sort(sorted_t* sd)
/*-------------------------------------------------------------*/
{
  free_mat(sd->mmat);
  free_mat_int(sd->midx);
  free(sd);
  return;
}

/*-------------------------------------------------------------*/
void copy_mat(mat_t *src, mat_t *dest)
/*-------------------------------------------------------------*/
{
  int i, j;
  for(i = 0; i < src->nr; ++i){
    for(j = 0; j < src->nc; ++j){
      dest->m[i][j] = src->m[i][j];
    }
  }
  return;
}

/*-------------------------------------------------------------*/
