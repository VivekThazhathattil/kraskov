#include <stdio.h>
#include "../data_handler.h"

int main(){
  mat_t* dat = get_data("../data/dset_test1_kraskov.h5", "/ds");
  int i, j;
  for(i = 0; i < 10; ++i){
    for(j = 0; j < 10; ++j){
      printf("%lf ", dat->m[i][j]);
    }
    printf("\n");
  }
  free_mat_t(dat);
  return 0;
}
