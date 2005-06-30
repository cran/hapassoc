/*  getWts.c   */

void getWts(int *size, int *ID, double *wts, double *num_prob) {
  int i, j;

  for (i = 0; i < size[0]; i++) {
    double sum_num_prob = 0;
    for (j = 0; j < size[0]; j++) {
      if (ID[j] == ID[i])
        sum_num_prob += num_prob[j];
    }
    wts[i] = num_prob[i]/sum_num_prob;
  }

}
