#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct MAT_s{
// nr: num rows, nc: num cols
  int nr, nc; 
  double **m;
} mat_t;

typedef struct SPEC_ELM_s{
  double val;
  int idx;
} spec_elm_t;

typedef struct MAT_INT_s{
// nr: num rows, nc: num cols
  int nr, nc; 
  int **m;
} mat_int_t;

typedef struct SORTED_s{
  struct MAT_s *mmat;
  struct MAT_INT_s *midx;
} sorted_t;

mat_t* init_mat(int nrows, int ncols){
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

mat_int_t* init_mat_int(int nrows, int ncols){
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

void print_mat(mat_t *mat){
  int i, j;
  for(i = 0; i < mat->nr; ++i){
    for(j = 0; j < mat->nc; ++j){
      printf("%lf ", mat->m[i][j]);
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

void free_mat(mat_t* mat){
  int i;
  for(i = 0; i < mat->nr; ++i){
    free(mat->m[i]);
  }
  free(mat->m);
  free(mat);
  return;
}

void free_mat_int(mat_int_t* mat){
  int i;
  for(i = 0; i < mat->nr; ++i){
    free(mat->m[i]);
  }
  free(mat->m);
  free(mat);
  return;
}

void free_mat_sort(sorted_t* sd){
  free_mat(sd->mmat);
  free_mat_int(sd->midx);
  free(sd);
  return;
}


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

void print_array(double *arr, int n){
  int i;
  for(i = 0; i < n; ++i){
    printf("%lf ", arr[i]);
  }
  printf("\n");
  return;
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

void copy_mat(mat_t *src, mat_t *dest){
  int i, j;
  for(i = 0; i < src->nr; ++i){
    for(j = 0; j < src->nc; ++j){
      dest->m[i][j] = src->m[i][j];
    }
  }
  return;
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

//lTODO:
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

int main(){
  double  x[3] = {1, 2, 3}; // 1st random variable 
  double  y[3] = {1, 5, 4}; // 2nd random variable
  int k = 3; // num nearest neighbors
  int n = 3; // num time snapshots

  remove_mean(x, n);
  remove_mean(y, n);

  std_normalize(x, n);
  std_normalize(y, n);

  printf("x:\n");
  print_array(x, n);
  printf("y:\n");
  print_array(y, n);

  // Get pdist between each time snapshot for a given time series

  double *pd_x = pdist(x, n); 
  double *pd_y = pdist(y, n);

  int n_pd = (n*(n-1))/2;
  printf("pd_x:\n");
  print_array(pd_x, n_pd);
  printf("pd_y:\n");
  print_array(pd_y, n_pd);

  double *pd = max_bw_two_array_elms(pd_x, pd_y, n_pd);
  printf("pd:\n");
  print_array(pd, n_pd);

  mat_t *pd_mat_x = square_form(pd_x, n_pd);
  mat_t *pd_mat_y = square_form(pd_y, n_pd);
  mat_t *pd_mat = square_form(pd, n_pd);

  printf("\n");
  printf("pd_mat_x:\n");
  print_mat(pd_mat_x);
  printf("pd_mat_y:\n");
  print_mat(pd_mat_y);
  printf("pd_mat:\n");
  print_mat(pd_mat);

  sorted_t *pd_mat_sort_x = sort_each_row(pd_mat_x);
  sorted_t *pd_mat_sort_y = sort_each_row(pd_mat_y);
  sorted_t *pd_mat_sort = sort_each_row(pd_mat);

  printf("pd_mat_sort_x->mmat:\n");
  print_mat(pd_mat_sort_x->mmat);
  printf("pd_mat_sort_x->midx:\n");
  print_mat_int(pd_mat_sort_x->midx);

  printf("pd_mat_sort_y->mmat:\n");
  print_mat(pd_mat_sort_y->mmat);
  printf("pd_mat_sort_y->midx:\n");
  print_mat_int(pd_mat_sort_y->midx);

  printf("pd_mat_sort->mmat:\n");
  print_mat(pd_mat_sort->mmat);
  printf("pd_mat_sort->midx:\n");
  print_mat_int(pd_mat_sort->midx);

  // free all the dynamically allocated variables
  free(pd_x);
  free(pd_y);
  free(pd);

  free_mat(pd_mat_x);
  free_mat(pd_mat_y);
  free_mat(pd_mat);

  free_mat_sort(pd_mat_sort_x);
  free_mat_sort(pd_mat_sort_y);
  free_mat_sort(pd_mat_sort);
  return 0;
}
