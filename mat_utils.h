#ifndef _MAT_UTILS_H
#define _MAT_UTILS_H

#include <stdio.h>
#include <stdlib.h>

typedef struct MAT_s{
	int nr, nc;
	double **m;	
} mat_t;

typedef struct MAT_INT_s{
// nr: num rows, nc: num cols
  int nr, nc; 
  int **m;
} mat_int_t;

typedef struct SORTED_s{
  struct MAT_s *mmat;
  struct MAT_INT_s *midx;
} sorted_t;

mat_t* init_mat(int, int);
mat_int_t* init_mat_int(int, int);
void free_mat(mat_t*);
void free_mat_int(mat_int_t*);
void free_mat_sort(sorted_t*);
void copy_mat(mat_t*, mat_t*);

#endif
