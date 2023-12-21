#ifndef _KRASKOV_H
#define _PRINT_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "oth_utils.h"
#include "mat_utils.h"
#include "data_handler.h"
#include "digamma.h"
#include "print_utils.h"

typedef struct SPEC_ELM_s{
  double val;
  int idx;
} spec_elm_t;

double mean(double* , int);
void remove_mean(double*, int);
double std(double*, int);
void std_normalize(double*, int);

double* pdist(double*, int);
double* max_bw_two_array_elms(double*, double*, int);
mat_t* square_form(double*, int);
void copy_mat_val_idx(double*, int, spec_elm_t*);
void create_index_mat(mat_int_t*);
int comp_fn1(const void*, const void*);
void spec_copy(spec_elm_t*, sorted_t*, int);
sorted_t* sort_each_row(mat_t*);
sorted_t* get_dist_mat_sort_xy(mat_t*, sorted_t*);
double* get_kth_dist(int k, sorted_t*);
void modify_eps(double*, double*, sorted_t*, sorted_t*, int);
int* get_dist_count(mat_t*, double*);
double* get_psi(int*, int);
double kraskov_mi(double*, double*, int, int, char);

#endif
