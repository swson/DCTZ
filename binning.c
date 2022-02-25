/**
 * @file binning.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief helper routines for binning mechanisms
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include "dctz.h"

#if 0
/* <= item < */
t_bin_id which_bin(double *bin_maxes, double item, double err_bound)
{
  /*
  double bin_width_side, bin_width_middle, whole_range;
  int index_middle = NBINS/2;
  int index_side = (NBINS - index_middle)/2;
  
  whole_range = NBINS*err_bound*2;
  bin_width_middle = (whole_range/3)/index_middle;
  bin_width_side = (whole_range/3)/index_side;
  */
#ifndef USE_BINARY_SEARCH
  int i;
#endif
  t_bin_id bin_id = NBINS; /* 0-{254,1022}: within range, 255/1023: out of range */

#ifdef USE_BINARY_SEARCH
  int first = 0, last = 254;
  int middle = (first+last)/2;

  while (first <= last) {
    if (item > bin_maxes[middle])
      first = middle + 1;
    else if (item < bin_maxes[middle] && item >= (bin_maxes[middle]-err_bound*2*BRSF))
      return (t_bin_id)middle;
    else 
      last = middle-1;
    middle = (first+last)/2;
  }
#else
  for (i=0; i<NBINS; i++)
    if (item < bin_maxes[i] && item >= (bin_maxes[i]-err_bound*2*BRSF))
      return (t_bin_id)i;
  /*
  for (i=0; i<index_side; i++) {
    if (item < bin_maxes[i] && item >= (bin_maxes[i]-bin_width_side*BRSF))
      return (t_bin_id)i;
  }
  for (i=index_side; i<index_side+index_middle; i++) {
    if (item < bin_maxes[i] && item >= (bin_maxes[i]-bin_width_middle*BRSF))
      return (t_bin_id)i;
  }
  for (i=index_side+index_middle; i<index_side+index_middle+index_side; i++){
    if (item < bin_maxes[i] && item >= (bin_maxes[i]-bin_width_side*BRSF))
      return (t_bin_id)i;
  }
*/
#endif

  return bin_id; /* not found */
}
#endif /* #if 0 */

void gen_bins(double min, double max, double *bin_maxes, double *bin_center, int nbins, double error_bound)
{
  double bin_width, range_min, range_max; 
  int i, half = nbins/2;
  /*
  double bin_width_side, bin_width_middle, whole_range;
  int index_middle = nbins/2;
  int index_side = (nbins - index_middle)/2;

  whole_range = nbins*error_bound*2;
  bin_width_middle = (whole_range/3)/index_middle;
  bin_width_side = (whole_range/3)/index_side;
  */
  bin_width = error_bound*2*BRSF;
  range_min = -(half*2+1)*(error_bound*BRSF);
  range_max = (half*2+1)*(error_bound*BRSF); /* FIX: set but not used */

#ifdef DEBUG
  printf("bin range: min = %e, max = %e\n", range_min, range_max);
#endif

  for (i=0; i<nbins; i++) {
    bin_maxes[i] = range_min + (i+1)*bin_width;
    bin_center[i] = bin_maxes[i] - (error_bound*BRSF);
  }
  /*
  for (i=0; i<index_side; i++){
    bin_maxes[i] = range_min + (i+1)*bin_width_side;
    bin_center[i] = bin_maxes[i] - (bin_width_side*BRSF)/2;
  }
  for (i=index_side; i<index_side+index_middle; i++){
    bin_maxes[i] = bin_maxes[index_side-1] + (i-index_side+1)*bin_width_middle;
    bin_center[i] = bin_maxes[i] - (bin_width_middle*BRSF)/2;
  }
  for (i=index_side+index_middle; i<index_side+index_middle+index_side; i++){
    bin_maxes[i] = bin_maxes[index_side+index_middle-1] + (i-index_side-index_middle+1)*bin_width_side;
    bin_center[i] = bin_maxes[i] - (bin_width_side*BRSF)/2;
  }
  */

  bin_center[half]=0.0; /* TODO: set to 0 for the center bin; otherwise, the computation above can't generate 0 for the center bin */

#ifdef DEBUG
  for (i=0; i<nbins; i++) 
    printf("bin_maxes[%d]: %e, ", i, bin_maxes[i]);
  printf("\n");

  for (i=0; i<nbins; i++) 
    printf("bin_center[%d]: %e, ", i, bin_center[i]);
  printf("\n");
#endif
}
