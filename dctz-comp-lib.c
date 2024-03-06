/**
 * @file dctz-comp-lib.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief DCTZ compression library routine
 * (C) 2019 University of Massachusetts Lowell.
       See LICENSE in the top-level directory.
*/

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#ifdef TIME_DEBUG
#include <sys/time.h>
#endif /* TIME_DEBUG */
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "zlib.h"
#include "dctz.h"
#include "dct.h"

#define DEF_MEM_LEVEL 8 // 8 (default), 9 uses maximum memory for optimal speed

t_bin_id conv_tbl[NBINS] = {
  254, 252, 250, 248, 246, 244, 242, 240, 238, 236, 234, 232, 230, 228, 226, 224,
  222, 220, 218, 216, 214, 212, 210, 208, 206, 204, 202, 200, 198, 196, 194, 192,
  190, 188, 186, 184, 182, 180, 178, 176, 174, 172, 170, 168, 166, 164, 162, 160,
  158, 156, 154, 152, 150, 148, 146, 144, 142, 140, 138, 136, 134, 132, 130, 128,
  126, 124, 122, 120, 118, 116, 114, 112, 110, 108, 106, 104, 102, 100, 98,  96,
  94,  92,  90,  88,  86,  84,  82,  80,  78,  76,  74,  72,  70,  68,  66,  64,
  62,  60,  58,  56,  54,  52,  50,  48,  46,  44,  42,  40,  38,  36,  34,  32,
  30,  28,  26,  24,  22,  20,  18,  16,  14,  12,  10,  8,   6,   4,   2,   0,
  1,   3,   5,   7,   9,   11,  13,  15,  17,  19,  21,  23,  25,  27,  29,  31,
  33,  35,  37,  39,  41,  43,  45,  47,  49,  51,  53,  55,  57,  59,  61,  63,
  65,  67,  69,  71,  73,  75,  77,  79,  81,  83,  85,  87,  89,  91,  93,  95,
  97,  99,  101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127,
  129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159,
  161, 163, 165, 167, 169, 171, 173, 175, 177, 179, 181, 183, 185, 187, 189, 191,
  193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223,
  225, 227, 229, 231, 233, 235, 237, 239, 241, 243, 245, 247, 249, 251, 253};

// 1, 2, 9, 17, 10, 3, 4, 11,
// 18, 25, 33, 26, 19, 12, 5, 6,
// 13, 20, 27, 34, 41, 49, 52, 35,
// 28, 21, 14, 7, 8, 15, 22, 29,
// 36, 43, 50, 57, 58, 51, 44, 37,
// 30, 23, 16, 24, 31, 38, 45, 52,
// 59, 60, 53, 46, 39, 32, 40, 47,
// 54, 61, 62, 55, 48, 56, 63, 64
t_bin_id zigzag[BLK_SZ] = {
  0,  1,  8,  16, 9,  2,  3,  10,
  17, 24, 32, 25, 18, 11, 4,  5,
  12, 19, 26, 33, 40, 48, 51, 34,
  27, 20, 13, 6,  7,  14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36,
  29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46,
  53, 60, 61, 54, 47, 55, 62, 63};

union 
{
  float *f;
  double *d;
} a, a_x;

union 
{
  float f;
  double d;
} bin_width, range_min, range_max, qt_factor;

void *compress_thread(void *arg)
{
  z_stream *defstream = (z_stream *)arg;
#ifdef DEBUG
  printf("compress started ...\n");
#endif
  deflate(defstream, Z_FINISH);
#ifdef DEBUG
  printf("done! compression...\n");
#endif
  uLong ret = defstream->total_out;
  deflateEnd(defstream);
  pthread_exit((void *)ret);
}

int dctz_compress(t_var *var, int N, size_t *outSize, t_var *var_z, double error_bound)
{
  int i, j, nblk, rem;
  //double TRUNC_BITS = (1.0/pow(error_bound, 1.55));
  //double TRUNC_BITS = pow(10.0, (ceil(abs(log2(error_bound)))/2.0)+2.0);
  //double TRUNC_BITS = pow(10.0, (ceil(abs(log2(error_bound))))-4.0);
  //float TRUNC_BITS = (1.0/(error_bound/0.01));
#ifdef TIME_DEBUG
  struct timeval start_t, end_t, gstart_t;
  double sf_t, dct_t, DC_AC_t, zlib_t, comp_t, malloc_t, genbin_t;
#endif
  t_bin_id *bin_index, *bin_indexz, *bin_indexz2;
#ifdef USE_TRUNCATE
  float *DC, *DCz, *DCz2, *AC_exact, *AC_exactz, *AC_exactz2;
#else
  double *DC, *DCz, *DCz2, *AC_exact, *AC_exactz, *AC_exactz2;
#endif
  struct header h;
  t_bstat bs;
  size_t type_size = 0; 
#ifdef USE_QTABLE
  union {
    double *d;
    float *f;
  } qtable;  /* Quantizer Table */
#endif

  if (var->datatype == DOUBLE)
    type_size = sizeof(double);
  else /* FLOAT */
    type_size = sizeof(float);

  if (var->datatype == DOUBLE) {
    if (NULL == (a_x.d = (double *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: a_x\n");
      exit(1);
    }
  }
  else { /* FLOAT */
    if (NULL == (a_x.f = (float *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: a_x\n");
      exit(1);
    }
  } 

  if (error_bound < 1E-6) {
    printf("ERROR BOUND is not acceptable");
    exit(1);
  }

#ifdef DEBUG
  for (i=0; i<BLK_SZ; i++) { /* show the first block */
    printf("a[%d] = %e\n", i, var->buf.d[i]); 
    if (i%BLK_SZ == 0 && i != 0) printf("\n");
  }
#endif

  if (NULL == (bin_index = (t_bin_id *)malloc(N*sizeof(t_bin_id)))) {
    fprintf(stderr, "Out of memory: bin_index[]\n");
    exit(1);
  }
  memset(bin_index, 0, sizeof(t_bin_id)*N);

#ifdef USE_QTABLE
  /* allocate space for qtable and initialize it */
  if (var->datatype == DOUBLE) {
    if (NULL == (qtable.d = (double *)malloc(BLK_SZ*sizeof(double)))) {
      fprintf(stderr, "Out of memory: qtable\n");
      exit(1);
    }
    for (i=0; i<BLK_SZ; i++) {
      qtable.d[i] = 0.0;
    }
  }
  else { /* FLOAT */
    if (NULL == (qtable.f = (float *)malloc(BLK_SZ*sizeof(float)))) {
      fprintf(stderr, "Out of memory: qtable\n");
      exit(1);
    }
    for (i=0; i<BLK_SZ; i++) {
      qtable.f[i] = 0.0;
    }
  }
#ifdef DEBUG
  for (i=0; i<BLK_SZ; i++) {
    printf("qtable[%d] = %e\n", i, (var->datatype == DOUBLE ? qtable.d[i] : qtable.f[i]));
  }
#endif /* DEBUG */
#endif /* USE_QTABLE */

#ifdef TIME_DEBUG
  gettimeofday(&start_t, NULL);
  gstart_t = start_t;
#endif
  
  /* determine scaling factor */
  calc_data_stat(var, &bs, N);

  if (var->datatype == DOUBLE) {
#ifdef DEBUG
    printf("scaling factor = %f\n", bs.sf.d);
#endif
    /* apply scaling factor */
    if (bs.sf.d != 1.0) {
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(bs.sf.d)
#endif
      for (i=0; i<N; i++) {
	var->buf.d[i] /= bs.sf.d;
	//var->buf.d[i] -= (bs.mean.d)/bs.sf.d;
      }
    }
  }
  else { /* FLOAT */
#ifdef DEBUG
    printf("scaling factor = %f\n", bs.sf.f);
#endif
    /* apply scaling factor */
    if (bs.sf.f != 1.0) {
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(bs.sf.f)
#endif
      for (i=0; i<N; i++) {
	var->buf.f[i] /= bs.sf.f;
	//var->buf.f[i] -= bs.mean.f/bs.sf.f;
      }
    }
  }
  
#ifdef TIME_DEBUG 
  gettimeofday(&end_t, NULL);
  sf_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));
  
  gettimeofday(&start_t, NULL);
#endif
  
  /* DCT over decomposed blocks */
  nblk = CEIL(N, BLK_SZ);
  rem = N % BLK_SZ;
#ifdef DEBUG
  printf("\nnumber of blocks = %d, remainder = %d\n", nblk, rem);
#endif

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

#ifdef USE_TRUNCATE
  if (NULL == (DCz = (float *)malloc(nblk*sizeof(float)))) {
    fprintf(stderr, "Out of memory: DCz[]\n");
    exit(1);
  }
  memset(DCz, 0, sizeof(float)*nblk); /* TODO: is it necessary? */
#else
  if (NULL == (DCz = (double *)malloc(nblk*sizeof(double)))) {
    fprintf(stderr, "Out of memory: DCz[]\n");
    exit(1);
  }
#endif

  if (NULL == (bin_indexz = (t_bin_id *)malloc(N*sizeof(t_bin_id)))) {
    fprintf(stderr, "Out of memory: bin_indexz[]\n");
    exit(1);
  }
  memset (bin_indexz, 0, sizeof(t_bin_id)*N);

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  malloc_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif /* TIME_DEBUG */

  int half = NBINS/2;
  if (var->datatype == DOUBLE) {
    bin_width.d = error_bound*2.0*BRSF;
    range_min.d = -(half*2+1)*(error_bound*BRSF);
    range_max.d = (half*2+1)*(error_bound*BRSF);
  }
  else { /* FLOAT */
    bin_width.f = error_bound*2.0*BRSF;
    range_min.f = -(half*2+1)*(error_bound*BRSF);
    range_max.f = (half*2+1)*(error_bound*BRSF);
  }

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  genbin_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

#ifdef USE_TRUNCATE
  if (NULL == (AC_exact = (float *)malloc(N*sizeof(float)))) {
    fprintf(stderr, "Out of memory: AC_exact\n");
    exit(1);
  }
  memset(AC_exact, 0, sizeof(float)*N); /* TODO: is it necessary? */
#else
  if (NULL == (AC_exact = (double *)malloc(N*sizeof(double)))) {
    fprintf(stderr, "Out of memory: AC_exact\n");
    exit(1);
  }
  memset(AC_exact, 0, sizeof(double)*N); /* TODO: is it necessary? */
#endif

#ifdef USE_TRUNCATE
  if (NULL == (AC_exactz = (float *)malloc(N*sizeof(float)))) {
    fprintf(stderr, "Out of memory: AC_exactz[]\n");
    exit(1);
  }
  memset(AC_exactz, 0, sizeof(float)*N); /* TODO: is it necessary? */
#else
  if (NULL == (AC_exactz = (double *)malloc(N*sizeof(double)))) {
    fprintf(stderr, "Out of memory: AC_exactz[]\n");
    exit(1);
  }
  memset(AC_exactz, 0, sizeof(double)*N); /* TODO: is it necessary? */
#endif

  if (var->datatype == DOUBLE)
    dct_init(BLK_SZ);
  else /* FLOAT */
    dct_init_f(BLK_SZ);

  int tot_AC_exact_count = 0;
  /* DCT block decomposed */
  for (i=0; i<nblk; i++) { /* for each decomposed blk */
    int l_blk_sz = ((i==nblk-1)&&(rem!=0))?rem:BLK_SZ;
    if ((i==nblk-1) && (rem!=0)) {
      if (var->datatype == DOUBLE) {
	dct_finish();
	dct_init(rem);
      }
      else { /* FLOAT */
	dct_finish_f();
	dct_init_f(rem);
      }
    }
    if (var->datatype == DOUBLE)
      dct_fftw(var->buf.d+i*BLK_SZ, a_x.d+i*BLK_SZ, l_blk_sz, nblk);
    else /* FLOAT */
      dct_fftw_f(var->buf.f+i*BLK_SZ, a_x.f+i*BLK_SZ, l_blk_sz, nblk);
    
#ifdef DEBUG
    printf("block %d: after DCT:\n", i);
    for (j=0; j<BLK_SZ && (i<3); j++){ /* show the first block only */
      printf("a_x[%d] = %e \n", i*BLK_SZ+j, a_x.d[i*BLK_SZ+j]);
    }
    printf("\n");
#endif
    
#ifdef USE_TRUNCATE
    DC[i] = (float)(var->datatype == DOUBLE ? a_x.d[i*BLK_SZ] : a_x.f[i*BLK_SZ]); /* save DC component in truncated */
#else
    DC[i] = (double)a_x.d[i*BLK_SZ]; /* save DC component */
#endif
#ifdef USE_QTABLE
    if (var->datatype == DOUBLE)
      qtable.d[0] = a_x.d[i*BLK_SZ];
    else /* FLOAT */
      qtable.f[0] = a_x.f[i*BLK_SZ];
#endif
    bin_index[i*BLK_SZ] = NBINS; /* store as it is */

    for (j=1; j<l_blk_sz; j++) {
      if (var->datatype == DOUBLE) {
	double item = a_x.d[i*BLK_SZ+j];
	t_bin_id bin_id;
	if (item < range_min.d || item > range_max.d) {
	  bin_id = NBINS;
#ifdef USE_QTABLE
	  /* get the highest values in each DCT location */
	  if (fabs(item) >= qtable.d[j])
	    qtable.d[j] = fabs(item);
#endif /* USE_QTABLE */
	}      
	else {
	  t_bin_id tmp_bin_id;
	  tmp_bin_id = (t_bin_id)((item-range_min.d)/bin_width.d);
	  bin_id = conv_tbl[tmp_bin_id];
	}
#ifdef DEBUG
	printf("bin_id = %d\n", bin_id);
#endif
	bin_index[i*BLK_SZ+j] = bin_id;
#ifdef DEBUG
	printf("a_x[%d]=%e => %d\n", i*BLK_SZ+j, (double)item, bin_id);
#endif
      } /* DOUBLE */
      else { /* FLOAT */
	//float item = fabsf(a_x.f[i*BLK_SZ+j]);
	float item = a_x.f[i*BLK_SZ+j];
	t_bin_id bin_id;
	if (item < (float)range_min.f || item > (float)range_max.f) {
	  bin_id = NBINS;
#ifdef USE_QTABLE
	  /* get the largest values in each DCT location */
	  if (fabsf(item) >= qtable.f[j])
	    qtable.f[j] = fabsf(item);
#endif /* USE_QTABLE */
	}      
	else {
	  t_bin_id tmp_bin_id;
	  tmp_bin_id = (t_bin_id)((item-range_min.f)/bin_width.f);
	  bin_id = conv_tbl[tmp_bin_id];
	}
#ifdef DEBUG
	printf("bin_id = %d\n", bin_id);
#endif /* DEBUG */
	bin_index[i*BLK_SZ+j] = bin_id;
#ifdef DEBUG
	printf("a_x[%d]=%e => %d\n", i*BLK_SZ+j, (float)item, bin_id);
#endif /* DEBUG */
	a_x.f[i*BLK_SZ+j] = item; /* update DCT coefficents with abs */
      } /* else: FLOAT */
    }
    /* end of of making quantization table  */
  }
  if (var->datatype == DOUBLE)
    dct_finish();
  else /* FLOAT */
    dct_finish_f();

#ifdef DCT_FILE_DEBUG
  FILE *fp = fopen("dct_result.bin", "w+");
  if (var->datatype == DOUBLE)
    fwrite(a_x.d, sizeof(double), N, fp);
  else /* FLOAT */
    fwrite(a_x.f, sizeof(float), N, fp);
  fclose(fp);

  fp = fopen("DC.bin", "w+");
  fwrite(DC, sizeof(float), nblk, fp);
  fclose(fp);
#endif /* DCT_FILE_DEBUG */

#ifdef USE_QTABLE
#ifdef DEBUG
  printf("Quantizer Table:\n");
  for (j=0; j<BLK_SZ ; j++) { 
    printf("before qtable[%d] = %e \n", j, (var->datatype == DOUBLE ? qtable.d[j] : qtable.f[j]));
  }
#endif /* DEBUG */

  FILE *fp_qtable = fopen("qtable.bin", "w+");
  if (var->datatype == DOUBLE)
    fwrite(qtable.d, sizeof(double), BLK_SZ, fp_qtable);
  else /* FLOAT */
    fwrite(qtable.f, sizeof(float), BLK_SZ, fp_qtable);
  fclose(fp_qtable);

  for (j=1; j<BLK_SZ; j++) {
    if (var->datatype == DOUBLE) {
      if (qtable.d[j] < 1.0) {
	qtable.d[j] = 1.0;
      }
    }
    else { /* FLOAT */
      if (qtable.f[j] < 1.0) {
	qtable.f[j] = 1.0;
      }
    }
  }
  
#ifdef DEBUG
  printf("Quantizer Table:\n");
  for (j=0; j<BLK_SZ ; j++) { 
    printf("after qtable[%d] = %e \n", j, (var->datatype == DOUBLE ? qtable.d[j] : qtable.f[j]));
  }
#endif /* DEBUG */
#endif /* USE_QTABLE */
  
#ifdef USE_QTABLE
  if (var->datatype == DOUBLE)
    qt_factor.d = (NBINS == 255 ? 10.0 : 2000.0);
  else /* FLAOT */
    qt_factor.f = (NBINS == 255 ? 10.0 : 2000.0);
#endif /* USE_QTABLE */
  
  for (i=0; i<nblk; i++) {
    int l_blk_sz = ((i==nblk-1)&&(rem != 0))?rem:BLK_SZ;
    for (j=1; j<l_blk_sz; j++) {
      t_bin_id bin_id;
      bin_id = bin_index[i*BLK_SZ+j];
      if (bin_id == NBINS) {
#ifdef USE_QTABLE
	if (var->datatype == DOUBLE) {
	  double item = a_x.d[i*BLK_SZ+j];
	  /* if out of bin area, normalize it to the area from range_max/range_min to range_max/range_min +/- error_bound */
	  if (item < range_min.d) {
	    item = (item/qtable.d[j])*error_bound*qt_factor.d + range_min.d;
	  } else if (item > range_max.d) {
	    item = (item/qtable.d[j])*error_bound*qt_factor.d + range_max.d;
	  }
	  a_x.d[i*BLK_SZ+j] = item; /* update a_x with updated value */
	  if (item < range_min.d || item > range_max.d) {
	    bin_id = NBINS;
#ifdef USE_TRUNCATE
	    AC_exact[tot_AC_exact_count++] = (float)a_x.d[i*BLK_SZ+j];
#else /* !USE_TRUNCATE */
	    AC_exact[tot_AC_exact_count++] = a_x.d[i*BLK_SZ+j];
#endif /* USE_TRUNCATE */
	  }
	  else {
	    t_bin_id tmp_bin_id;
	    tmp_bin_id = (t_bin_id)((item-range_min.d)/bin_width.d);
	    bin_id = conv_tbl[tmp_bin_id];
	  }
#ifdef DEBUG
	  printf("a_x[%d]=%e => %d\n", i*BLK_SZ+j, (double)item, bin_id);
#endif /* DEBUG */
	} /* if DOUBLE */
	else { /* FLOAT */
	  float item = a_x.f[i*BLK_SZ+j];
	  /* if out of bin area, normalize it to the area from range_max/range_min to range_max/range_min +/- error_bound */
	  if (item < range_min.f) {
	    item = (item/qtable.f[j])*error_bound*qt_factor.f + range_min.f;
	  } else if (item > range_max.f) {
	    item = (item/qtable.f[j])*error_bound*qt_factor.f + range_max.f;
	  }
	  a_x.f[i*BLK_SZ+j] = item; /* update a_x with updated value */
	  if (item < range_min.f || item > range_max.f) {
	    bin_id = NBINS;
	    //AC_exact[tot_AC_exact_count++] = roundf(a_x.f[i*BLK_SZ+j]*TRUNC_BITS)/TRUNC_BITS;
	    AC_exact[tot_AC_exact_count++] = a_x.f[i*BLK_SZ+j];
	  }
	  else {
	    t_bin_id tmp_bin_id;
	    tmp_bin_id = (t_bin_id)((item-range_min.f)/bin_width.f);
	    bin_id = conv_tbl[tmp_bin_id];
	  }
#ifdef DEBUG
	  printf("a_x[%d]=%e => %d\n", i*BLK_SZ+j, item, bin_id);
#endif /* DEBUG */
	}
#else /* !USE_QTABLE, i.e., EC */
#ifdef USE_TRUNCATE
	//AC_exact[tot_AC_exact_count++] = (var->datatype == DOUBLE ? (roundf((float)a_x.d[i*BLK_SZ+j])*TRUNC_BITS)/TRUNC_BITS : roundf(a_x.f[i*BLK_SZ+j]*TRUNC_BITS)/TRUNC_BITS);
	AC_exact[tot_AC_exact_count++] = (var->datatype == DOUBLE ? (float)a_x.d[i*BLK_SZ+j] : a_x.f[i*BLK_SZ+j]);
#else
	AC_exact[tot_AC_exact_count++] = (var->datatype == DOUBLE ? a_x.d[i*BLK_SZ+j] : a_x.f[i*BLK_SZ+j]);
#endif /* USE_TRUNCATE */
#endif /* !USE_QTABLE, i.e., EC */
    } /* bin_id == NBINS */
  } /* for j */
} /* for i */

#ifdef SIZE_DEBUG
  printf("total AC_exact_count = %d (%.2f percent) \n", tot_AC_exact_count, (double)(tot_AC_exact_count)/(double)N*100.0);
#endif

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  dct_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

#ifdef DEBUG
  int bin_freq[NBINS+1] = {0};

  i=0;
  while (i < N) {
    bin_freq[(int)bin_index[i++]]++;
  }
  printf("i=%d\n", i);
  
  int sum = 0;
  printf("bin_freq: ");
  for (i=0; i<NBINS+1; i++) {
    printf("%d, ", bin_freq[i]);
    sum += bin_freq[i];
  }
  printf("sum=%d\n", sum);
#endif

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);

  DC_AC_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

  char bin_index_file[640];
  FILE *fp_index;
  sprintf(bin_index_file, "bin_index.bin"); 
  fp_index = fopen(bin_index_file, "wb");
  fwrite(bin_index, N, 1, fp_index);
  fclose(fp_index);

  char AC_exact_file[640];
  FILE *fp_AC_exact;
  sprintf(AC_exact_file, "AC_exact.bin"); 
  fp_AC_exact = fopen(AC_exact_file, "wb");
  fwrite(AC_exact, tot_AC_exact_count*sizeof(float), 1, fp_AC_exact);
  fclose(fp_AC_exact);

#if 0
  /* zigzag scanning */
  for (i=0; i<nblk; i++)
    for (j=0; j<BLK_SZ; j++)
      bin_index[i*BLK_SZ+j] = zigzag[bin_index[i*BLK_SZ+j]];
#endif
  
#ifdef DEBUG
  printf("tot_AC_exact_count=%d\n", tot_AC_exact_count);
#ifdef USE_QTABLE  
  printf("bin_index before compression = %lu\n", N*sizeof(t_bin_id)); // TODO: k should be N? k was used before AC normalization (HPEC'20)
#else  
  printf("bin_index before compression = %lu\n", N*sizeof(t_bin_id));
#endif  
#ifdef USE_TRUNCATE
  printf("DC before compression = %lu\n", nblk*sizeof(float));
  printf("AC_exact before compression = %lu\n", tot_AC_exact_count*sizeof(float));
#else
  printf("DC before compression = %lu\n", nblk*sizeof(double));
  printf("AC_exact before compression = %lu\n", tot_AC_exact_count*sizeof(double));
#endif
#endif

  pthread_t thread[3];
  pthread_attr_t attr;            /* thread attributes (left at defaults) */

  /* set defaults (not all pthread implementations default to joinable) */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  /* setup for compress */
  z_stream defstream[3];
  
  defstream[0].zalloc = Z_NULL;
  defstream[0].zfree = Z_NULL;
  defstream[0].opaque = Z_NULL;
  
  /* compress bin_index */
#ifdef USE_QTABLE
  uLong ucompSize_binindex = N*sizeof(t_bin_id);
#else 
  uLong ucompSize_binindex = N*sizeof(t_bin_id);
#endif  
  uLong compSize_binindex = compressBound(ucompSize_binindex);

  int windowBits = 15; 
  deflateInit2(&defstream[0], Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);
 
  defstream[0].avail_in = ucompSize_binindex;
  defstream[0].next_in = (Bytef *)bin_index;
  defstream[0].avail_out = compSize_binindex;
  defstream[0].next_out = (Bytef *)bin_indexz;
  defstream[0].data_type = Z_UNKNOWN; /* Z_TEXT(Z_ASCII), Z_BINARY, Z_UNKNOWN */

  if (pthread_create(&thread[0], &attr, compress_thread, (void *)&defstream[0])) {
    fprintf(stderr, "Error creating thread\n");
    exit(0);
  }
  
  /* compress DC */
  defstream[1].zalloc = Z_NULL;
  defstream[1].zfree = Z_NULL;
  defstream[1].opaque = Z_NULL;

#ifdef USE_TRUNCATE
  uLong ucompSize_DC = nblk*sizeof(float);
  uLong compSize_DC = compressBound(ucompSize_DC);
#else
  uLong ucompSize_DC = nblk*sizeof(double);
  uLong compSize_DC = compressBound(ucompSize_DC);
#endif

  deflateInit2(&defstream[1], Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);
 
  defstream[1].avail_in = ucompSize_DC;
  defstream[1].next_in = (Bytef *)DC;
  defstream[1].avail_out = compSize_DC;
  defstream[1].next_out = (Bytef *)DCz;
  defstream[1].data_type = Z_UNKNOWN;

  if (pthread_create(&thread[1], &attr, compress_thread, (void *)&defstream[1])) {
    fprintf(stderr, "Error creating thread\n");
    exit(0);
  }

  /* compress AC_exact */
  defstream[2].zalloc = Z_NULL;
  defstream[2].zfree = Z_NULL;
  defstream[2].opaque = Z_NULL;
  
#ifdef USE_TRUNCATE
  uLong ucompSize_AC_exact = tot_AC_exact_count*sizeof(float);
  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#else
  uLong ucompSize_AC_exact = tot_AC_exact_count*sizeof(double);
  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#endif

  deflateInit2(&defstream[2], Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);

  defstream[2].avail_in = ucompSize_AC_exact;
  defstream[2].next_in = (Bytef *)AC_exact;
  defstream[2].avail_out = compSize_AC_exact;
  defstream[2].next_out = (Bytef *)AC_exactz;
  defstream[2].data_type = Z_UNKNOWN;

  if (pthread_create(&thread[2], &attr, compress_thread, (void *)&defstream[2])) {
    fprintf(stderr, "Error creating thread\n");
    exit(0);
  }

#ifdef USE_TRUNCATE
  //uLong ucompSize_AC_exact = tot_AC_exact_count*sizeof(float);
  //  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#else
  //uLong ucompSize_AC_exact = tot_AC_exact_count*sizeof(double);
  // uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#endif
  void *ret;
  for (i=0; i<3; i++) {
    pthread_join(thread[i], &ret);
#ifdef DEBUG
    printf("thread %d joined\n", i);
#endif
    switch (i) {
    case 0:
      compSize_binindex = (uLong)ret;
      break;
    case 1:
      compSize_DC = (uLong)ret;
      break;
    case 2:
      compSize_AC_exact = (uLong)ret;
      break;
    }
  }

  pthread_attr_destroy(&attr);
#if 0
  compSize_binindex = defstream[0].total_out; /* update with actual size */
  deflateEnd(&defstream[0]);

  compSize_DC = defstream[1].total_out; /* update with actual size */
  deflateEnd(&defstream[1]);

  compSize_AC_exact_count = defstream[2].total_out; /* update with actual size */
  deflateEnd(&defstream[2]);
#endif

  bin_indexz2 = (t_bin_id*)realloc(bin_indexz, compSize_binindex); /* TODO: check error */

#ifdef SIZE_DEBUG
  printf("Compressed bin_index size is: %lu\n", compSize_binindex);
#endif

  DCz2 = realloc(DCz, compSize_DC); /* TODO: check error */
#ifdef SIZE_DEBUG
  printf("Compressed DC size is: %lu\n", compSize_DC);
#endif

  AC_exactz2 = realloc(AC_exactz, compSize_AC_exact); /* TODO: check error */
#ifdef SIZE_DEBUG
  printf("Compressed AC_exact size is: %lu\n", compSize_AC_exact);
#endif

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  double comp_rate;

  zlib_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(start_t.tv_sec*1000000 + start_t.tv_usec));
  comp_t = (double)((end_t.tv_sec*1000000 + end_t.tv_usec)-(gstart_t.tv_sec*1000000 + gstart_t.tv_usec));
  comp_rate = (N*sizeof(double)/(double)(1024*1024))/(comp_t/1000000);

  printf("sf_t=%f(s), dct_t=%f(s), zlib_t(compress)=%f(s)\n", sf_t/1000000, dct_t/1000000, zlib_t/1000000);
  printf("malloc_t=%f(s), genbin=%f(s), DC_AC_t=%f(s)\n", malloc_t/1000000, genbin_t/1000000, DC_AC_t/1000000);
  printf("comp_time = %f (s), compression rate = %f (MB/s)\n", comp_t/1000000, comp_rate);
#endif

  *outSize = sizeof(struct header) + compSize_binindex + compSize_DC + compSize_AC_exact;
#ifdef USE_QTABLE
  *outSize += (var->datatype == DOUBLE ? sizeof(double)*BLK_SZ : sizeof(float)*BLK_SZ);
#endif

  h.datatype = var->datatype;
  h.num_elements = N;
  h.error_bound = error_bound;
  h.tot_AC_exact_count = tot_AC_exact_count;
  if (var->datatype == DOUBLE) {
    h.scaling_factor.d = bs.sf.d;
    h.mean.d = bs.mean.d;
  }
  else { /* FLOAT */ 
    h.scaling_factor.f = bs.sf.f;
    h.mean.f = bs.mean.f;
  }
  h.bindex_sz_compressed = compSize_binindex;
  h.DC_sz_compressed = compSize_DC;
  h.AC_exact_sz_compressed = compSize_AC_exact;
#ifdef USE_QTABLE
  h.bindex_count = N;
#endif  
  //h.AC_exact_count_sz_compressed = compSize_AC_exact_count;

  unsigned char *cur_p;
  if (var->datatype == DOUBLE)
    cur_p = (unsigned char *)(var_z->buf.d);
  else
    cur_p = (unsigned char *)(var_z->buf.f);
  memcpy(cur_p, &h, sizeof(struct header));
  cur_p += sizeof(struct header);
  memcpy(cur_p, bin_indexz2, compSize_binindex);
  cur_p += compSize_binindex;
  //memcpy (cur_p, AC_exact_countz2, compSize_AC_exact_count);
  //cur_p += compSize_AC_exact_count;
  memcpy(cur_p, DCz2, compSize_DC);
  cur_p += compSize_DC;
  memcpy(cur_p, AC_exactz2, compSize_AC_exact);
  cur_p += compSize_AC_exact;
#ifdef USE_QTABLE
  if (var->datatype == DOUBLE)
    memcpy(cur_p, qtable.d, BLK_SZ*sizeof(double));
  else /* FLAOT */
    memcpy(cur_p, qtable.f, BLK_SZ*sizeof(float));
#endif /* USE_QTABLE */

  if (var->datatype == DOUBLE) {
    free(a_x.d);
#ifdef USE_QTABLE
    free(qtable.d);
#endif
  }
  else { /* FLOAT */
    free(a_x.f);
#ifdef USE_QTABLE
    free(qtable.f);
#endif
  }
  free(DC);
  free(DCz2); 
  free(AC_exact);
  free(AC_exactz2);
  free(bin_index);
  free(bin_indexz2);
  
#ifndef SIZE_DEBUG
  printf("outSize = %zu\n", *outSize);
#endif

  return(1);
}
