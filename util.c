/**
 * @file util.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief utility routines
 * (C) 2019 University of Massachusetts Lowell.
       See LICENSE in the top-level directory.
*/

#include "dctz.h"

void calc_data_stat(t_var *in, t_bstat *bs, int N)
{
  int i;

  if (in->datatype == DOUBLE) {
    double sum = 0.0;
    bs->max.d = fabs(in->buf.d[0]);
    bs->min.d = fabs(in->buf.d[0]);

    for (i=1; i<N; i++) {
      if (fabs(in->buf.d[i]) > bs->max.d) bs->max.d = fabs(in->buf.d[i]);
      if (fabs(in->buf.d[i]) < bs->min.d) bs->min.d = fabs(in->buf.d[i]);
      sum += in->buf.d[i];
    }

    bs->mean.d = sum/N;
    bs->sf.d = pow(10, ceil(log10(bs->max.d)) - SF_ADJ_AMT);
  }
  else { /* FLOAT */
    float sum = 0.0;
    bs->max.f = fabsf(in->buf.f[0]);
    bs->min.f = fabsf(in->buf.f[0]);

    for (i=1; i<N; i++) {
      if (fabsf(in->buf.f[i]) > bs->max.f) bs->max.f = fabsf(in->buf.f[i]);
      if (fabsf(in->buf.f[i]) < bs->min.f) bs->min.f = fabsf(in->buf.f[i]);
      sum += in->buf.f[i];
    }

    bs->mean.f = sum/N;
    bs->sf.f = powf(10, ceil(log10f(bs->max.f)) - SF_ADJ_AMT);
  }
}

int iceil(double d)
{
  int i = (int)d;
  if (d == (double)i)
    return i;
  return i + 1;
}

double calc_psnr(t_var *var, t_var *var_r, int N, double error_bound)
{
  double sum_sq=0.0, mse, rmse, relative_range, psnr, min, max, relative_error;
  double maxdiff = 0.0;
  int i;

  if (var->datatype == DOUBLE) {
    min = max = var->buf.d[0];

    for (i=1; i<N; i++) {
      if (var->buf.d[i] > max) max = var->buf.d[i];
      if (var->buf.d[i] < min) min = var->buf.d[i];
    }
    for (i=0; i<N; i++) {
      if(fabs(var->buf.d[i]-var_r->buf.d[i])>maxdiff){
        maxdiff = fabs(var->buf.d[i]-var_r->buf.d[i]);
      }
      double error = var->buf.d[i]-var_r->buf.d[i];
      sum_sq += (error*error);
    }
  }
  else { /* FLOAT */
    min = max = var->buf.f[0];

    for (i=1; i<N; i++) {
      if (var->buf.f[i] > max) max = var->buf.f[i];
      if (var->buf.f[i] < min) min = var->buf.f[i];
    }
    for (i=0; i<N; i++) {
      if(fabs(var->buf.f[i]-var_r->buf.f[i])>maxdiff){
        maxdiff = fabs(var->buf.f[i]-var_r->buf.f[i]);
      }
      float error = var->buf.f[i]-var_r->buf.f[i];
      sum_sq += (error*error);
    }
  }
  mse = sum_sq/N;
  rmse = sqrt(mse);
  relative_range = max - min;
  psnr = 20*log10(relative_range/rmse);
  relative_error = maxdiff/relative_range;
  printf("Max relative error = %.6f\n", relative_error);
#if 0
  if (relative_error < error_bound){
    printf("Guarantee the error bound\n");
  }else{
    printf("Not guarantee the error bound\n");
  }
#endif
  return psnr;
}
