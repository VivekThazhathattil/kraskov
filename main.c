/* Ref: https://github.com/stefgrs/Mutual-Information-script/blob/master/MIpair2.m */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_utils.h"
#include "data_handler.h"
#include "digamma.h"
#include "print_utils.h"

typedef struct SPEC_ELM_s{
  double val;
  int idx;
} spec_elm_t;

double vabs(double num){
  if(num > 0){
    return num;
  }
  return -num;
}

double vmax(double n1, double n2){
  if(n1 > n2){
    return n1;
  }
  return n2;
}

double mean(double *arr, int n){
  int i;
  double sum = 0.0;
  for(i = 0; i < n; ++i){
    sum += arr[i];
  }
  return sum/n;
}

void remove_mean(double *arr, int n){
  double m = mean(arr, n);
  int i;
  for(i = 0; i < n; ++i){
    arr[i] = arr[i] - m;
  }
  return;
}

double std(double *arr, int n){
  int i;
  double m = mean(arr, n);
  double sum = 0;
  for(i = 0; i < n; ++i){
    sum += (arr[i] - m)*(arr[i] - m);
  }
  if(n == 1){
    printf("ERROR: Cannot normalize when n = 1!\n");
    return -1;
  }
  return sqrt(sum/(n - 1));
}

void std_normalize(double *arr, int n){
  double d = std(arr, n);
  int i;
  for(i = 0; i < n; ++i){
    arr[i] = arr[i]/d;
  }
}

double* pdist(double *arr, int n){
  int i, j, k = 0;
  int m = ((n - 1)*n)/2;
  double *res = (double*) calloc(m, sizeof(double));
  for(i = 0; i < n; ++i){
    for(j = i + 1; j < n; ++j){
      res[k] = vabs(arr[i] - arr[j]);
      ++k;
    }
  }
  return res;
}

double* max_bw_two_array_elms(double *a, double *b, int n){
  int i;
  double* c = (double*) malloc(sizeof(double) * n);
  for(i = 0; i < n; ++i){
    c[i] = vmax(a[i], b[i]);
  }
  return c;
}

mat_t* square_form(double* arr, int n){
  mat_t* mat = init_mat(n, n);
  int i, j, k = 0;
  for(i = 0; i < mat->nr; ++i){
    for(j = i; j < mat->nc; ++j){
      if(i == j){
        mat->m[i][j] = 0.0;
      }
      else{
        mat->m[i][j] = arr[k];
        ++k;
      }
    }
  }
  for(i = 0; i < mat->nr; ++i){
    for(j = 0; j < i; ++j){
      mat->m[i][j] = mat->m[j][i];
    }
  }
  return mat;
}

void copy_mat_val_idx(double *src, int n, spec_elm_t* elms){
  int i;
  for(i = 0; i < n; ++i){
    elms[i].val = src[i];
    elms[i].idx = i;
  }
  return;
}

void create_index_mat(mat_int_t *mat){
  int i, j;
  for(i = 0; i < mat->nr; ++i){
    for(j = 0; j < mat->nc; ++j){
      mat->m[i][j] = j;
    }
  }
  return;
}

int comp_fn1(const void* l, const void* r){
  double res = ((spec_elm_t*)l)->val - ((spec_elm_t*)r)->val;
  if(res < 0){
    return -1;
  }
  else if(res == 0){
    return 0;
  }
  return 1;
}

void spec_copy(spec_elm_t* elms, sorted_t* sd, int idx){
  int j; 
  for(j = 0; j < (sd->mmat)->nc; ++j){
    //printf("%lf %d\n", elms[j].val, elms[j].idx);
    (sd->mmat)->m[idx][j] = elms[j].val;
    (sd->midx)->m[idx][j] = elms[j].idx;
  }
  return;
}

sorted_t* sort_each_row(mat_t* mat){
  int i;
  sorted_t* sd = (sorted_t*) malloc(sizeof(sorted_t));
  spec_elm_t* elms = 
    (spec_elm_t*) malloc(sizeof(spec_elm_t) * mat->nc);
  sd->mmat = init_mat(mat->nr, mat->nc);
  sd->midx = init_mat_int(mat->nr, mat->nc);
  for(i = 0; i < mat->nr; ++i){
    copy_mat_val_idx(mat->m[i], mat->nc, elms);
    qsort(elms, mat->nc, sizeof(spec_elm_t), comp_fn1);
    spec_copy(elms, sd, i);
  }
  free(elms);
  return sd;
}

sorted_t* get_dist_mat_sort_xy(mat_t* mat, sorted_t* sd){
  sorted_t* res = (sorted_t*) malloc(sizeof(sorted_t));
  res->mmat = init_mat(mat->nr, mat->nc);
  res->midx = init_mat_int(mat->nr, mat->nc);
  mat_int_t *tm = (sd->midx);
  int i, j;
  for(i = 0; i < tm->nr; ++i){
    for(j = 0; j < tm->nc; ++j){
      (res->mmat)->m[i][j] = mat->m[i][tm->m[i][j]];
      (res->midx)->m[i][j] = tm->m[i][j];
    }
  }
  return res;
}

double* get_kth_dist(int k, sorted_t* sd){
  if(k >= (sd->mmat)->nr){
    printf("ERROR: Number of time snapshots cannot be <= k \n");
    exit(1);
  }
  double* res = 
    (double*) malloc(sizeof(double) * (sd->mmat)->nr);
  int i;
  for(i = 0; i < (sd->mmat)->nr; ++i){
    res[i] = (sd->mmat)->m[i][k];
  }
  return res;
}

void modify_eps(double* eps_x, double* eps_y, 
    sorted_t* sd_x, sorted_t* sd_y, int k){
 int n = (sd_x->mmat)->nr, i, j;
 for(i = 0; i < n; ++i){
  int flag_y = 1;
  for(j = k - 1; j >= 1; --j){
    if((sd_x->mmat)->m[i][j] > eps_x[i]){
      eps_x[i] = (sd_x->mmat)->m[i][j];
      flag_y = 0;
    }
  }
  if(flag_y){
    for(j = k - 1; j >= 1; --j){
      if((sd_y->mmat)->m[i][j] > eps_y[i]){
        eps_y[i] = (sd_y->mmat)->m[i][j];
      }
    }
  }
 }
 return;
}

int* get_dist_count(mat_t* dmat, double* eps){
  int* count = (int*) malloc(sizeof(int) * dmat->nr);
  int i, j, sum;
  for(i = 0; i < dmat->nr; ++i){
    sum = 0;
    for(j = 0; j < dmat->nc; ++j){
      if(dmat->m[i][j] <= eps[i]){
        ++sum;
      }
    }
    count[i] = sum - 1;
    if(count[i] < 0){
      count[i] = 0;
    }
  }
  return count;
}

double* get_psi(int* vals, int n){
  int i;
  double* res = (double*) malloc(sizeof(double) * n);
  for(i = 0; i < n; ++i){
    res[i] = (double) psi((double)vals[i]);
  }
  return res;
}

int main(){
  /* h5ls should return output of dimensions `nsamples x nfeatures` */
  mat_t* dat = get_data("data/dset_test1_kraskov.h5", "/ds");
  double* x = dat->m[0];
  double* y = dat->m[1];

  int k = 3; // num nearest neighbors
  int n = dat->nc; // num time snapshots

  printf("no. rows = %d, no. cols = %d\n", dat->nr, dat->nc);

  remove_mean(x, n);
  remove_mean(y, n);

  std_normalize(x, n);
  std_normalize(y, n);

  //printf("x:\n");
  //print_array(x, n);
  //printf("y:\n");
  //print_array(y, n);

  // Get pdist between each time snapshot for a given time series

  double *pd_x = pdist(x, n); 
  double *pd_y = pdist(y, n);

  int n_pd = (n*(n-1))/2;
  //printf("pd_x:\n");
  //print_array(pd_x, n_pd);
  //printf("pd_y:\n");
  //print_array(pd_y, n_pd);

  double *pd = max_bw_two_array_elms(pd_x, pd_y, n_pd);
  //printf("pd:\n");
  //print_array(pd, n_pd);
  
  printf("[VIVEK]2a\n");

  mat_t *dist_mat_x = square_form(pd_x, n);
  mat_t *dist_mat_y = square_form(pd_y, n);
  mat_t *dist_mat = square_form(pd, n);

  printf("[VIVEK]3\n");
  //printf("\n");
  //printf("dist_mat_x:\n");
  //print_mat(dist_mat_x);
  //printf("dist_mat_y:\n");
  //print_mat(dist_mat_y);
  //printf("dist_mat:\n");
  //print_mat(dist_mat);

  sorted_t *dist_mat_sort = sort_each_row(dist_mat);
  sorted_t *dist_mat_sort_x = 
    get_dist_mat_sort_xy(dist_mat_x, dist_mat_sort);
  sorted_t *dist_mat_sort_y = 
    get_dist_mat_sort_xy(dist_mat_y, dist_mat_sort);

  //printf("dist_mat_sort_x->mmat:\n");
  //print_mat(dist_mat_sort_x->mmat);
  //printf("dist_mat_sort_x->midx:\n");
  //print_mat_int(dist_mat_sort_x->midx);
  printf("[VIVEK]4\n");

  //printf("dist_mat_sort_y->mmat:\n");
  //print_mat(dist_mat_sort_y->mmat);
  //printf("dist_mat_sort_y->midx:\n");
  //print_mat_int(dist_mat_sort_y->midx);

  //printf("dist_mat_sort->mmat:\n");
  //print_mat(dist_mat_sort->mmat);
  //printf("dist_mat_sort->midx:\n");
  //print_mat_int(dist_mat_sort->midx);

  double *eps_x = get_kth_dist(k, dist_mat_sort_x);
  double *eps_y = get_kth_dist(k, dist_mat_sort_y);

  printf("[VIVEK]5\n");
  //printf("eps_x:\n");
  //print_array(eps_x, n);
  //printf("eps_y:\n");
  //print_array(eps_y, n);

  modify_eps(eps_x, eps_y, dist_mat_sort_x, 
      dist_mat_sort_y, k);

  //printf("eps_x (modified):\n");
  //print_array(eps_x, n);
  //printf("eps_y (modified):\n");
  //print_array(eps_y, n);
  printf("[VIVEK]6\n");

  int* nx = get_dist_count(dist_mat_x, eps_x);
  int* ny = get_dist_count(dist_mat_y, eps_y);

  //printf("nx:\n");
  //print_array_int(nx, n);
  //printf("ny:\n");
  //print_array_int(ny, n);
  printf("[VIVEK]7\n");

  double* psi_vals_x = get_psi(nx, n);
  double* psi_vals_y = get_psi(ny, n);
  double mean_psi_x = mean(psi_vals_x, n);
  double mean_psi_y = mean(psi_vals_y, n);

  printf("[VIVEK]8\n");
  printf("mean_psi_x: %0.4lf, mean_psi_y: %0.4lf\n", mean_psi_x, mean_psi_y);

  double mi = psi(k) - (1.0/(double)k) - mean_psi_x - mean_psi_y + psi(n);
  printf("MI = %lf\n", mi);

  // free all the dynamically allocated variables
  free(pd_x);
  free(pd_y);
  free(pd);

  free_mat(dist_mat_x);
  free_mat(dist_mat_y);
  free_mat(dist_mat);

  free_mat_sort(dist_mat_sort_x);
  free_mat_sort(dist_mat_sort_y);
  free_mat_sort(dist_mat_sort);

  free(eps_x);
  free(eps_y);

  free(nx);
  free(ny);

  free(psi_vals_x);
  free(psi_vals_y);

  return 0;
}
