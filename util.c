/**
 * @file util.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief utility routines
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include "dctz.h"

void calc_data_stat (double *in, struct bstat *bs, int N)
{
  int i;

  bs->max = fabs (in[0]);
  bs->min = fabs (in[0]);

  for (i=1; i<N; i++) {
    if (fabs (in[i]) > bs->max) bs->max = fabs (in[i]);
    if (fabs (in[i]) < bs->min) bs->min = fabs (in[i]);
  }

  bs->sf = ceil(log10(bs->max));
}

void calc_data_stat_f (float *in, struct bstat *bs, int N)
{
  int i;

  bs->max = fabs (in[0]);
  bs->min = fabs (in[0]);

  for (i=1; i<N; i++) {
    if (fabs (in[i]) > bs->max) bs->max = fabs (in[i]);
    if (fabs (in[i]) < bs->min) bs->min = fabs (in[i]);
  }
}
