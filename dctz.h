/**
 * @file dctz.h
 * @author Seung Woo Son
 * @date July 2019
 * @brief header for DCTZ routines
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#ifndef _DCTZ_H_
#define _DCTZ_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fftw3.h"
#include "zlib.h"
#include "dct.h"

#define DCTZ_VERSION "0.2.1"
#define DCTZ_VERSION_MAJOR 0
#define DCTZ_VERSION_MINOR 2
#define DCTZ_VERSION_PATCH 1

#define BLK_SZ 64
#define BRSF 1 /* bin range scaling factor: 1: no scaling */

#define MAX(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

#define MIN(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a < _b ? _a : _b; })

#define CEIL(a,b) (a+b-1)/b

typedef enum {
  FLOAT = 0,
  DOUBLE
} t_datatype;

typedef struct
{
  t_datatype datatype;
  double err_bound;
  char *var_name;
  union
  {
    float *f;
    double *d;
  } buf;
} t_var;

/* unsigned char => 1 byte: 0-255 */
/* unsigned short => 2 byte: 0-1023 */
typedef unsigned char t_bin_id;

#define NBITS (sizeof(t_bin_id)<<3) /* # of bits (8 or 16) for representing bin index */
#define NBINS ((1 << (NBITS)) - 1)

typedef struct {
  union
  {
    double d;
    float f;
  } min;
  union
  {
    double d;
    float f;
  } max;
  union
  {
    double d;
    float f;
  } range;
  union
  {
    double d;
    float f;
  } sf;
} t_bstat;

struct header 
{
  t_datatype datatype;
  unsigned int num_elements;
  double error_bound;
  unsigned int tot_AC_exact_count;
  union
  {
    double d;
    float f;
  } scaling_factor;
  unsigned int bindex_sz_compressed;
  unsigned int DC_sz_compressed;
  unsigned int AC_exact_sz_compressed;
#ifdef USE_QTABLE 
  unsigned int bindex_count;
#endif  
  /*  unsigned int AC_exact_count_sz_compressed; */
};

void calc_data_stat(t_var *in, t_bstat *bs, int N);
int ceili(double i);
void gen_bins(double min, double max, double *bin_maxes, double *bin_center, int nbins, double error_bound);
void *compress_thread(void *arg);
int dctz_compress(t_var *var, int N, size_t *outSize, t_var *var_z, double error_bound);
int dctz_decompress(t_var *var_z, t_var *var_r);
double calc_psnr(t_var *var, t_var *var_r, int N);
#endif
