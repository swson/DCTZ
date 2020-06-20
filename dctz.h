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

#define DCTZ_VERSION "0.2.0"
#define DCTZ_VERSION_MAJOR 0
#define DCTZ_VERSION_MINOR 2
#define DCTZ_VERSION_PATCH 0

#define BLK_SZ 64
#define NBITS 16 /* # of bits (8 or 16) for representing bin index */
#define NBINS ((1 << (NBITS)) - 1)
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

enum dtype {data_type_double, data_type_float};

union buffers {
  float *bf;
  double *bd;
} bt;

union {
  unsigned char onebyte;
  unsigned short twobyte;
} Bin_Id;

struct bstat {
  double min;
  double max;
  double range;
  double sf;
};

struct header 
{
  unsigned int num_elements;
  double error_bound;
  unsigned int tot_AC_exact_count;
  double scaling_factor;
  unsigned int bindex_sz_compressed;
  unsigned int DC_sz_compressed;
  unsigned int AC_exact_sz_compressed;
#ifdef USE_QTABLE 
  unsigned int bindex_count;
#endif  
//  unsigned int AC_exact_count_sz_compressed;
};

void calc_data_stat (double *in, struct bstat *bs, int N);
void calc_data_stat_f (float *in, struct bstat *bs, int N);
int ceili (double i);
void gen_bins (double min, double max, double *bin_maxes, double *bin_center, int nbins, double error_bound);
void *compress_thread (void *arg);
int dctz_compress (double *d, int N, size_t *outSize, char *a_z, double error_bound);
int dctz_compress_float (float *f, int N, size_t *outSize, char *a_z, double error_bound);
int dctz_decompress (char *, double *);
int dctz_decompress_float (char *, float *);

#endif
