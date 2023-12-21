#include "print_utils.h"

void print_mat(mat_t *mat){
  int i, j;
  for(i = 0; i < mat->nr; ++i){
    for(j = 0; j < mat->nc; ++j){
      printf("%0.4lf ", mat->m[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  return; 
}

void print_mat_int(mat_int_t *mat){
  int i, j;
  for(i = 0; i < mat->nr; ++i){
    for(j = 0; j < mat->nc; ++j){
      printf("%d ", mat->m[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  return; 
}

void print_array(double *arr, int n){
  int i;
  for(i = 0; i < n; ++i){
    printf("%0.4lf ", arr[i]);
  }
  printf("\n");
  return;
}

void print_array_int(int *arr, int n){
  int i;
  for(i = 0; i < n; ++i){
    printf("%d ", arr[i]);
  }
  printf("\n");
  return;
}
