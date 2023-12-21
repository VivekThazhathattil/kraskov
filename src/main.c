/* Ref: https://github.com/stefgrs/Mutual-Information-script/blob/master/MIpair2.m */

#include <stdio.h>
#include "kraskov.h"

/*-------------------------------------------------------------*/
int main()
/*-------------------------------------------------------------*/
{
  /* h5ls should return output of dimensions `nsamples x nfeatures` */
  mat_t* dat0 = get_data("data/dset_test2_kraskov.h5", "/ds");
  mat_t* dat1 = get_data("data/dset_test2_kraskov_0.h5", "/ds");
  double* y = dat1->m[0];

  int i, k = 3; // num nearest neighbors
  int n = dat0->nc; // num time snapshots
                    
  printf("Rows = %d, Cols = %d \n", dat0->nr, dat0->nc);

  for(i = 0; i < dat0->nr; ++i){
    double* x = dat0->m[i];
    kraskov_mi(x, y, n, k);
  }

  free_mat_t(dat0);
  free_mat_t(dat1);

  return 0;
}
