#include "oth_utils.h"

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
