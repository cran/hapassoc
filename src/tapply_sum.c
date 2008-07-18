/*  tapply_sum.c   */

void tapply_sum (int *size, int *factor, double *vector, int *ID, double *sum) {
  int i, j;
  
  for (i = 0; i < *size; i++) {
     j = 0;
     while (ID[i] != factor[j]) j++;
     sum[j] += vector[i];
  }
}
