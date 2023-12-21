/* Ref: https://github.com/stefgrs/Mutual-Information-script/blob/master/MIpair2.m */

#include <omp.h>
#include <stdio.h>
#include "kraskov.h"

#define PRINT_FACTOR 1000

/*-------------------------------------------------------------*/
int main()
/*-------------------------------------------------------------*/
{
  /* h5ls should return output of dimensions `nsamples x nfeatures` */

  /* Source points */
  mat_t* dat0 = get_data("data/dset_test2_kraskov.h5", "/ds");

  /* Target points */
  mat_t* dat1 = get_data("data/dset_test2_kraskov_0.h5", "/ds");

  double* tgt = dat1->m[0];
  int i, k = 3; // num nearest neighbors
  int n = dat0->nc; // num time snapshots 
  int n_features = dat0->nr; // num features
  int num_threads = 24;
                    
  printf("Rows = %d, Cols = %d \n", dat0->nr, dat0->nc);

  double* vals = (double*) malloc(sizeof(double) * n_features);
  for(i = 0; i < n_features; ++i){
    vals[i] = 0;
  }

  int counter = 0;

  omp_set_num_threads(num_threads);

  /* normalize the values beforehand if running in parallel */
  for(i = 0; i < n_features; ++i){
    remove_mean(dat0->m[i], n);
    std_normalize(dat0->m[i], n);
  }
  remove_mean(tgt, n);
  std_normalize(tgt, n);

  #pragma omp parallel for firstprivate(counter)
  for(i = 0; i < dat0->nr; ++i){
    ++counter;
    double* src = dat0->m[i];
    vals[i] = kraskov_mi(src, tgt, n, k, 'n');
    if(counter % PRINT_FACTOR == 0){
      printf("[%d/~%d]\n", 
          counter, n_features/num_threads);
    }
  }

  save_vec_to_h5(vals, n_features, "data/out/test2_result.h5");

  free_mat_t(dat0);
  free_mat_t(dat1);
  free(vals);

  return 0;
}

/*-------------------------------------------------------------*/
