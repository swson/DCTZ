/**
 * @file binning.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief helper routines for binning mechanisms
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include "dctz.h"

/* <= item < */
unsigned char which_bin (double *bin_maxes, double item, double err_bound)
{
#ifndef USE_BINARY_SEARCH
  int i;
#endif
  unsigned char bin_id = 255; /* 0-254: within range, 255: out of range */

#ifdef USE_BINARY_SEARCH
  int first = 0, last = 254;
  int middle = (first+last)/2;

  while (first <= last) {
    if (item > bin_maxes[middle])
      first = middle + 1;
    else if (item < bin_maxes[middle] && item >= (bin_maxes[middle]-err_bound*2*BRSF))
      return (unsigned char)middle;
    else 
      last = middle-1;
    middle = (first+last)/2;
  }
#else
  for (i=0; i<NBINS; i++)
    if (item < bin_maxes[i] && item >= (bin_maxes[i]-err_bound*2*BRSF))
      return (unsigned char)i;
#endif

  return bin_id; /* not found */
}

void gen_bins (double min, double max, double *bin_maxes, double *bin_center, int nbins, double error_bound)
{
  double bin_width, range_min, range_max; 
  int i, half = nbins/2;

  bin_width = error_bound*2*BRSF;
  range_min = -(half*2+1)*(error_bound*BRSF);
  range_max = (half*2+1)*(error_bound*BRSF); // FIX: set but not used

#ifdef DEBUG
  printf ("bin range: min = %e, max = %e\n", range_min, range_max);
#endif

  for (i=0; i<nbins; i++) {
    bin_maxes[i] = range_min + (i+1)*bin_width;
    bin_center[i] = bin_maxes[i] - (error_bound*BRSF);
  }
  bin_center[half]=0.0; /* TODO: set to 0 for the center bin; otherwise, the computation above can't generate 0 for the center bin */

#ifdef DEBUG
  for (i=0; i<nbins; i++) 
    printf ("bin_maxes[%d]: %e, ", i, bin_maxes[i]);
  printf ("\n");

  for (i=0; i<nbins; i++) 
    printf ("bin_center[%d]: %e, ", i, bin_center[i]);
  printf ("\n");
#endif
}
