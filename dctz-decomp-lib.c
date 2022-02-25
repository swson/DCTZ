/**
 * @file dctz-decomp-lib.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief DCTZ decompression library routine
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <math.h>
#ifdef TIME_DEBUG
#include <sys/time.h> /* struct timeval */
#endif
#include "fftw3.h"
#include "zlib.h"
#include "dctz.h"
#include "dct.h"

union
{
  float *f;
  double *d;
} var_xr; /* var to store the result of IDCT */

int dctz_decompress(t_var *var_z, t_var *var_r)
{
  int N, i, j, nblk, rem;
#ifdef TIME_DEBUG
  struct timeval start_t, end_t, gstart_t; /* for measuring timing */
  double sf_t, zlib_t, idct_t, decomp_t, malloc_t, genbin_t;
#endif
  double *bin_maxes, *bin_center;
  t_bin_id *bin_index, *bin_indexz;
  unsigned int tot_AC_exact_count;
#ifdef USE_TRUNCATE
  float *DC, *DCz, *AC_exact, *AC_exactz;
#else
  double *DC, *DCz, *AC_exact, *AC_exactz;
#endif /* USE_TRUNCATE */
  struct header h;
  double error_bound;
#ifdef USE_QTABLE
  double *qtable;    /* Quantizer Table */
#endif /* USE_QTABLE */
  t_datatype datatype = var_z->datatype;

  unsigned char *cur_p;
  if (datatype == DOUBLE) {
    cur_p = (unsigned char *)var_z->buf.d;
    memcpy(&h, var_z->buf.d, sizeof(struct header));
    cur_p = (unsigned char *)var_z->buf.d + sizeof(struct header);
  }
  else { /* FLOAT */
    cur_p = (unsigned char *)var_z->buf.f;
    memcpy(&h, var_z->buf.f, sizeof(struct header));
    cur_p = (unsigned char *)var_z->buf.f + sizeof(struct header);
  }
  
  N = h.num_elements;
  error_bound = h.error_bound;
  nblk = CEIL(N, BLK_SZ);
  rem = N % BLK_SZ;
  tot_AC_exact_count = h.tot_AC_exact_count;
#ifdef USE_QTABLE
  int K = h.bindex_count;
  if (NULL == (bin_indexz = (t_bin_id *)malloc(K*sizeof(t_bin_id)))) {
    fprintf(stderr, "Out of memory: bin_indexz[]\n");
    exit(1);
  }
#else
  if (NULL == (bin_indexz = (t_bin_id *)malloc(N*sizeof(t_bin_id)))) {
    fprintf(stderr, "Out of memory: bin_indexz[]\n");
    exit(1);
  }

#endif  

#ifdef DEBUG
  printf("nitems=%d, tot_AC_exact_count=%d, scaling_factor=%e, bindex_sz_compressed=%d, DC_sz_compressed=%d, AC_exact_sz_compressed=%d,\n", h.num_elements, h.tot_AC_exact_count, (datatype == DOUBLE)?h.scaling_factor.d:h.scaling_factor.f, h.bindex_sz_compressed, h.DC_sz_compressed, h.AC_exact_sz_compressed);
#endif
#ifdef DEBUG
  printf("N=%d, nblk=%d, SF=%e\n", N, nblk, (datatype == DOUBLE) ? h.scaling_factor.d : h.scaling_factor.f);
#endif

#ifdef USE_TRUNCATE
  if (NULL == (DCz = (float *)malloc(nblk*sizeof(float)))) {
    fprintf(stderr, "Out of memory: DCz[]\n");
    exit(1);
  }
#else
  if (NULL == (DCz = (double *)malloc(nblk*sizeof(double)))) {
    fprintf(stderr, "Out of memory: DCz[]\n");
    exit(1);
  }
#endif

#ifdef USE_TRUNCATE
  if (NULL == (AC_exactz = (float *)malloc(N*sizeof(float)))) {
    fprintf(stderr, "Out of memory: AC_exactz[]\n");
    exit(1);
  }
#else
  if (NULL == (AC_exactz = (double *)malloc(N*sizeof(double)))) {
    fprintf(stderr, "Out of memory: AC_exactz[]\n");
    exit(1);
  }
#endif

#ifdef USE_TRUNCATE
  if (NULL == (AC_exact = (float *)malloc(tot_AC_exact_count*sizeof(float)))) {
    fprintf(stderr, "Out of memory: AC_exact\n");
    exit(1);
  }
#else
  if (NULL == (AC_exact = (double *)malloc(tot_AC_exact_count*sizeof(double)))) {
    fprintf(stderr, "Out of memory: AC_exact[]\n");
    exit(1);
  }
#endif

#ifdef USE_QTABLE
  /* Initialize Quantizer Table */
  if (NULL == (qtable = (double *)malloc(BLK_SZ*sizeof(double)))) {
      fprintf(stderr, "Out of memory: qtable\n");
      exit(1);
   }
#endif

  if (datatype == DOUBLE) {
    if (NULL == (var_xr.d = (double *)malloc(N*sizeof(double)))) {
      fprintf(stderr, "Out of memory: var_xr\n");
      exit(1);
    }
  }
  else { /* FLOAT */
    if (NULL == (var_xr.f = (float *)malloc(N*sizeof(float)))) {
      fprintf(stderr, "Out of memory: var_xr\n");
      exit(1);
    }
  }

  memcpy(bin_indexz, cur_p, h.bindex_sz_compressed);
  cur_p += h.bindex_sz_compressed;
//  memcpy(AC_exact_countz, cur_p, h.AC_exact_count_sz_compressed);
//  cur_p += h.AC_exact_count_sz_compressed;
  memcpy(DCz, cur_p, h.DC_sz_compressed);
  cur_p += h.DC_sz_compressed;
  memcpy(AC_exactz, cur_p, h.AC_exact_sz_compressed);
#ifdef USE_QTABLE
  cur_p += h.AC_exact_sz_compressed;
  memcpy(qtable, cur_p, BLK_SZ*sizeof(double));
#endif

#ifdef USE_QTABLE
#ifdef DEBUG
    printf("Quantizer Table:\n");
    for (j=1; j<BLK_SZ ; j++) { /* Show Quantizer Table */
      printf("qtable[%d] = %e \n", j, qtable[j]);
    }
#endif
#endif

#ifdef TIME_DEBUG
  gettimeofday(&start_t, NULL);
  gstart_t = start_t;
#endif

#ifdef USE_QTABLE
  uLong ucompSize_binindex = K*sizeof(t_bin_id);
  if (NULL == (bin_index = (t_bin_id *)malloc(K*sizeof(t_bin_id)))) {
    fprintf(stderr, "Out of memory: bin_index[]\n");
    exit(1);
  }
#else  
  uLong ucompSize_binindex = N*sizeof(t_bin_id);
  if (NULL == (bin_index = (t_bin_id *)malloc(N*sizeof(t_bin_id)))) {
    fprintf(stderr, "Out of memory: bin_index[]\n");
    exit(1);
  }
#endif  
  uLong compSize_binindex = compressBound(ucompSize_binindex);

  uLong ucompSize_DC = nblk*sizeof(float); /*  float or double? */
  uLong compSize_DC = compressBound(ucompSize_DC);

  uLong ucompSize_AC_exact = N*sizeof(float); /* tot_ac_exact instead of N? */
  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);

  /* setup for decompress */
  z_stream infstream[3];
  infstream[0].zalloc = Z_NULL;
  infstream[0].zfree = Z_NULL;
  infstream[0].opaque = Z_NULL;

  inflateInit(&infstream[0]);

  /* decompress bin_index */
  infstream[0].avail_in = compSize_binindex;
  infstream[0].next_in = (Bytef *)bin_indexz;
  infstream[0].avail_out = ucompSize_binindex;
  infstream[0].next_out = (Bytef *)bin_index;

  inflate(&infstream[0], Z_NO_FLUSH);
  inflateEnd(&infstream[0]);
#ifndef SIZE_DEBUG
  printf("uncompressed bin_index size is: %lu\n", infstream[0].total_out);
#endif

  /* decompress DC */
  infstream[1].zalloc = Z_NULL;
  infstream[1].zfree = Z_NULL;
  infstream[1].opaque = Z_NULL;

#ifdef USE_TRUNCATE
  if (NULL == (DC = (float *)malloc(nblk*sizeof(float)))) {
    fprintf(stderr, "Out of memory: DC[]\n");
    exit(1);
  }
#else
  if (NULL == (DC = (double *)malloc(nblk*sizeof(double)))) {
    fprintf(stderr, "Out of memory: DC[]\n");
    exit(1);
  }
#endif

  inflateInit(&infstream[1]);

  infstream[1].avail_in = compSize_DC;
  infstream[1].next_in = (Bytef *)DCz;
  infstream[1].avail_out = ucompSize_DC;
  infstream[1].next_out = (Bytef *)DC;

  inflate(&infstream[1], Z_NO_FLUSH);
  inflateEnd(&infstream[1]);
#ifndef SIZE_DEBUG
  printf("uncompressed DC size is: %lu\n", infstream[1].total_out);
#endif

  /* decompress AC_exact */
  infstream[2].zalloc = Z_NULL;
  infstream[2].zfree = Z_NULL;
  infstream[2].opaque = Z_NULL;
  
#ifdef USE_TRUNCATE
  if (NULL == (AC_exact = (float *)malloc(N*sizeof(float)))) {
    fprintf(stderr, "Out of memory: AC_exact[]\n");
    exit(1);
  }
#else
  if (NULL == (AC_exact = (double *)malloc(N*sizeof(double)))) {
    fprintf(stderr, "Out of memory: AC_exact[]\n");
    exit(1);
  }
#endif

  inflateInit(&infstream[2]);

  infstream[2].avail_in = compSize_AC_exact;
  infstream[2].next_in = (Bytef *)AC_exactz;
  infstream[2].avail_out = ucompSize_AC_exact;
  infstream[2].next_out = (Bytef *)AC_exact;

  inflate(&infstream[2], Z_NO_FLUSH);
  inflateEnd(&infstream[2]);
#ifndef DEBUG
  printf("uncompressed AC_exact size is: %lu\n", infstream[2].total_out);
#endif

  free(bin_indexz);
  free(DCz);
  free(AC_exactz);

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  zlib_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif
  
  /* restore AC_exact */
  nblk = CEIL(h.num_elements, BLK_SZ);
 
  if (NULL == (bin_maxes = (double *)malloc(NBINS*sizeof(double)))) {
    fprintf(stderr, "Out of memory: bin_maxes\n");
    exit(1);
  }

  if (NULL == (bin_center = (double *)malloc(NBINS*sizeof(double)))) {
    fprintf(stderr, "Out of memory: bin_center\n");
    exit(1);
  }

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  malloc_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif
  
  gen_bins(0, 0, bin_maxes, bin_center, NBINS, error_bound);

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  genbin_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif
  
  unsigned int pos = 0;
  unsigned int c = N;

  double range_max = error_bound * NBINS;
  double range_min = -error_bound * NBINS;
  double qt_factor = (NBINS == 255 ? 10.0 : 2000.0);

  if (datatype == DOUBLE)
    dct_init(BLK_SZ);
  else /* FLOAT */
    dct_init_f(BLK_SZ);
  
  /* IDCT block decomposed */
  for (i=0; i<nblk; i++) { /* for each decomposed blk */
    int l_blk_sz = ((i==nblk-1)&&(rem != 0))?rem:BLK_SZ;
    if (datatype == DOUBLE) {
      var_xr.d[i*BLK_SZ] = DC[i]; /* if USE_TRUNCATE, then float -> double */
#ifdef DEBUG
      printf("var_xr[%d]=%e\n", i*BLK_SZ, var_xr.d[i*BLK_SZ]);
#endif
      for (j=1; j<l_blk_sz; j++) {
#ifdef USE_QTABLE
	t_bin_id sbin_id;
#endif
	if (bin_index[i*BLK_SZ+j] == NBINS) {
#ifdef USE_QTABLE /* DCT_QT */
	  sbin_id = bin_index[c++];
	  
	  if (sbin_id == NBINS) {
	    var_xr.d[i*BLK_SZ+j] = AC_exact[pos]; /* if USE_TRUNCATE, then float -> double */
	    pos++;
	  }
	  else {
	    var_xr.d[i*BLK_SZ+j] = bin_center[sbin_id];
	  }
	  if (var_xr.d[i*BLK_SZ+j] > 0) {
	    var_xr.d[i*BLK_SZ+j] = ((var_xr.d[i*BLK_SZ+j] - range_max)/(error_bound*qt_factor))*qtable[j];
	  } else {
	    var_xr.d[i*BLK_SZ+j] = ((var_xr.d[i*BLK_SZ+j] - range_min)/(error_bound*qt_factor))*qtable[j];
	  }
#else /* DCT_EC */
	  var_xr.d[i*BLK_SZ+j] = AC_exact[pos]; /* if USE_TRUNCATE, then float -> double */
	  pos++;
#endif /* USE_QTABLE */
	}
	else {
	  var_xr.d[i*BLK_SZ+j] = bin_center[bin_index[i*BLK_SZ+j]];
	}
#ifdef DEBUG
	printf("after var_xr[%d]=%e\n", i*BLK_SZ+j, var_xr.d[i*BLK_SZ]+j);
#endif
      }

      if ((i==nblk-1) && (rem!=0)) {
	dct_finish();
	dct_init(rem);
      }
      
      ifft_idct(l_blk_sz, var_xr.d+i*BLK_SZ, var_r->buf.d+i*BLK_SZ);

#ifdef DEBUG
      printf("block %d: after IDCT:\n", i);
      for (j=0; j<BLK_SZ && (i < 3); j++) { /* show the first block results */
	printf("var_r[%d] = %e \n", i*BLK_SZ+j, var_r->buf.d[i*BLK_SZ+j]);
      }
#endif
    } /* DOUBLE */
    else { /* FLOAT */
      var_xr.f[i*BLK_SZ] = DC[i]; /* if USE_TRUNCATE, then float -> double */
#ifdef DEBUG
      printf("var_xr[%d]=%e\n", i*BLK_SZ, var_xr.f[i*BLK_SZ]);
#endif
      for (j=1; j<l_blk_sz; j++) {
#ifdef USE_QTABLE
	t_bin_id sbin_id;
#endif
	if (bin_index[i*BLK_SZ+j] == NBINS) {
#ifdef USE_QTABLE /* DCT_QT */
	  sbin_id = bin_index[c++];
	  
	  if (sbin_id == NBINS) {
	    var_xr.f[i*BLK_SZ+j] = AC_exact[pos]; /* if USE_TRUNCATE, then float -> double */
	    pos++;
	  }
	  else {
	    var_xr.f[i*BLK_SZ+j] = bin_center[sbin_id];
	  }
	  if (var_xr.d[i*BLK_SZ+j] > 0) {
	    var_xr.f[i*BLK_SZ+j] = ((var_xr.f[i*BLK_SZ+j] - range_max)/(error_bound*qt_factor))*qtable[j];
	  } else {
	    var_xr.f[i*BLK_SZ+j] = ((var_xr.f[i*BLK_SZ+j] - range_min)/(error_bound*qt_factor))*qtable[j];
	  }
#else /* DCT_EC */
	  var_xr.f[i*BLK_SZ+j] = AC_exact[pos]; /* if USE_TRUNCATE, then float -> double */
	  pos++;
#endif /* USE_QTABLE */
	}
	else {
	  var_xr.f[i*BLK_SZ+j] = bin_center[bin_index[i*BLK_SZ+j]];
	}
#ifdef DEBUG
	printf("after var_xr[%d]=%e\n", i*BLK_SZ+j, var_xr.f[i*BLK_SZ]+j);
#endif
      }

      if ((i==nblk-1) && (rem!=0)) {
	dct_finish_f();
	dct_init_f(rem);
      }

      ifft_idct_f(l_blk_sz, var_xr.f+i*BLK_SZ, var_r->buf.f+i*BLK_SZ);

#ifdef DEBUG
      printf("block %d: after IDCT:\n", i);
      for (j=0; j<BLK_SZ && (i < 3); j++) { /* show the first block results */
	printf("var_r[%d] = %e \n", i*BLK_SZ+j, var_r->buf.f[i*BLK_SZ+j]);
      }
#endif
    } /* FLOAT */
  }
  
  free(DC);

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  idct_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));
  
  gettimeofday(&start_t, NULL);
#endif

  if (datatype == DOUBLE) {
    double xscale = pow(10, (h.scaling_factor.d)-1);
    /*  deapply scaling factor */
    if (h.scaling_factor.d != 1.0) {
      for (i=0; i<N; i++)
	var_r->buf.d[i] *= xscale;
    }
  }
  else { /* FLOAT */
    float xscale = pow(10, (h.scaling_factor.f)-1);
    /* deapply scaling factor */
    if (h.scaling_factor.f != 1.0) {
      for (i=0; i<N; i++)
	var_r->buf.f[i] *= xscale;
    }
  }

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  sf_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&end_t, NULL);

  double decomp_rate;

  zlib_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));
  decomp_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(gstart_t.tv_sec*1000000 + gstart_t.tv_usec));
  decomp_rate = (N*sizeof(double)/(double)(1024*1024))/(decomp_t/1000000);

  printf("sf_t=%f(s), idct_t=%f(s), zlib_t(uncompress)=%f(s)\n", sf_t/1000000, idct_t/1000000, zlib_t/1000000);
  printf("malloc_t=%f(s), genbin=%f(s)\n", malloc_t/1000000, genbin_t/1000000);
  printf("decomp_time = %f (s), decompression rate = %f (MB/s)\n", decomp_t/1000000, decomp_rate);
#endif

  if (datatype == DOUBLE)
    idct_finish();
  else /* FLOAT */
    idct_finish_f();
  free(AC_exact);
  free(bin_index);
  free(bin_maxes);
  free(bin_center);
#ifdef USE_QTABLE
  free(qtable);
#endif /* USE_QTABLE */

  return(1);
}
