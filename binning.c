/**
 * @file binning.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief helper routines for binning mechanisms
 * (C) 2019 University of Massachusetts Lowell.
       See LICENSE in the top-level directory.
*/

#include "dctz.h"

void gen_bins(double min, double max, double *bin_center, int nbins, double error_bound)
{
  double bin_width;
  int i;

  bin_width = error_bound*2*BRSF;

  bin_center[0] = 0.0;
  for (i=1; i<nbins; i++) {
    int tmp_i = (i%2) ? ((i/2)+1) : -(i/2); 
    bin_center[i] = tmp_i*bin_width;
  }

#ifdef DEBUG
  for (i=0; i<nbins; i++) 
    printf("bin_center[%d]: %e, ", i, bin_center[i]);
  printf("\n");
#endif
}

void gen_bins_f(float min, float max, float *bin_center, int nbins, float error_bound)
{
  float bin_width;
  int i;

  bin_width = error_bound*2*BRSF;

  bin_center[0] = 0.0;
  for (i=1; i<nbins; i++) {
    int tmp_i = (i%2) ? ((i/2)+1) : -(i/2); 
    bin_center[i] = tmp_i*bin_width;
  }

#ifdef DEBUG
  for (i=0; i<nbins; i++) 
    printf("bin_center[%d]: %e, ", i, bin_center[i]);
  printf("\n");
#endif
}
