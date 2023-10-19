/**
 * @file dct-test.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief testing DCT performance
 * (C) 2019 University of Massachusetts Lowell.
       See LICENSE in the top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include "dctz.h"

/* gcc -o dct-test dct-test.c dct.c dct-float.c util.c -lfftw3 -lfftw3f -lm -Wall -g */

int main(int argc, char * argv[])
{
  size_t r4=0,r3=0,r2=0,r1=0;
  size_t typesize = 0;
  char *oriFilePath, outputFilePath[640];
  double *d, *d_x, *d_r;
  float *f, *f_x, *f_r;
  int datatype;
  int N, i, nblk, rem;
  
  if (argc < 4) {
    printf("Test case: %s -d|-f [srcFilePath] [dimension sizes...]\n", argv[0]);
    printf("Example: %s -d testdata/x86/testfloat_8_8_128.dat 8 8 128 dctz-ec(1E-3) \n", argv[0]);

    exit(0);
  }
  
  assert(argc >= 4);

  if (argc >= 4) { /* 1D */
    r1 = N = atoi(argv[3]);
  }
  if (argc >= 5) { /* 2D */
    r2 = atoi(argv[4]);
    N = r1*r2;
  }
  if (argc >= 6) { /* 3D */
    r3 = atoi(argv[5]);
    N = r1*r2*r3;
  } 
  if (argc >= 7) { /* 4D */
    r4 = atoi(argv[6]);
    N = r1*r2*r3*r4;
  }
  
  printf("total number = %d\n", N);

  nblk = CEIL(N, BLK_SZ);
  rem = N % BLK_SZ;
  printf("nblk=%d, rem=%d\n", nblk, rem);

  oriFilePath = argv[2];
  
  FILE *fp_in = fopen(oriFilePath, "rb");
  if (fp_in == NULL) {
    perror("Failed: ");
    printf("File Not Found\n");
    return(1);
  }
  
  if (!strcmp(argv[1], "-d")) { /* double */
    typesize = sizeof(double);
    datatype = DOUBLE;
    if (NULL == (d = (double *)malloc(N*typesize))) {
      fprintf(stderr, "Out of memory: a\n");
      exit(1);
    }

    if (NULL == (d_x = (double *)malloc(N*typesize))) {
      fprintf(stderr, "Out of memory: a_z\n");
      exit(1);
    }
    fread(d, typesize, N, fp_in);
    dct_init(BLK_SZ);
    for (i=0; i<nblk; i++) {
      int l_blk_sz = ((i==nblk-1)&&(rem!=0))?rem:BLK_SZ;
      if ((i==nblk-1)&&(rem!=0)) {
	dct_finish();
	dct_init(rem);
      }
      dct_fftw(d+i*BLK_SZ, d_x+i*BLK_SZ, l_blk_sz, nblk);
    }
    dct_finish();
  } 
  else { /* float */
    typesize = sizeof(float);
    datatype = FLOAT;
    if (NULL == (f = (float *)malloc(N*typesize))) {
      fprintf(stderr, "Out of memory: a\n");
      exit(1);
    }
    if (NULL == (f_x = (float *)malloc(N*typesize))) {
      fprintf(stderr, "Out of memory: a\n");
      exit(1);
    }
    fread(f, typesize, N, fp_in);
    dct_init_f(BLK_SZ);
    for (i=0; i<nblk; i++) {
      int l_blk_sz = ((i==nblk-1)&&(rem!=0))?rem:BLK_SZ;
      if ((i==nblk-1)&&(rem!=0)) {
	dct_finish_f();
	dct_init_f(rem);
      }
      dct_fftw_f(f+i*BLK_SZ, f_x+i*BLK_SZ, l_blk_sz, nblk);
    }
    dct_finish_f();
  }	  
  
  fclose(fp_in);

  FILE *fp_x;
  int icount;
  size_t outSize = N*typesize; /* outSize is same as input */
  
  sprintf(outputFilePath, "%s.x", oriFilePath);

  fp_x = fopen(outputFilePath, "w");
  if (datatype == DOUBLE)
    icount = fwrite(d_x, outSize, 1, fp_x);
  else /* float */
    icount = fwrite(f_x, outSize, 1, fp_x);
  if (icount != 1) {
    printf("Write DCT coefficients to a file failed: %lu != %d!\n", outSize, icount);
    exit(1);
  }
  fclose(fp_x);

  /* IDCT */
  if (!strcmp(argv[1], "-d")) { /* double */
    typesize = sizeof(double);
    datatype = DOUBLE;

    if (NULL == (d_r = (double *)malloc(N*typesize))) {
      fprintf(stderr, "Out of memory: a_z\n");
      exit(1);
    }
    dct_init(BLK_SZ);
    for (i=0; i<nblk; i++) {
      int l_blk_sz = ((i==nblk-1)&&(rem!=0))?rem:BLK_SZ;
      if ((i==nblk-1)&&(rem!=0)) {
	dct_finish();
	dct_init(rem);
      }
      ifft_idct(l_blk_sz, d_x+i*BLK_SZ, d_r+i*BLK_SZ);
    }
    dct_finish();
  } 
  else { /* float */
    typesize = sizeof(float);
    datatype = FLOAT;

    if (NULL == (f_r = (float *)malloc(N*typesize))) {
      fprintf(stderr, "Out of memory: a\n");
      exit(1);
    }
    dct_init_f(BLK_SZ);
    for (i=0; i<nblk; i++) {
      int l_blk_sz = ((i==nblk-1)&&(rem!=0))?rem:BLK_SZ;
      if ((i==nblk-1)&&(rem!=0)) {
	dct_finish_f();
	dct_init_f(rem);
      }
      ifft_idct_f(l_blk_sz, f_x+i*BLK_SZ, f_r+i*BLK_SZ);
    }
    dct_finish_f();
  }	  

  FILE *fp_r;
  sprintf(outputFilePath, "%s.r", oriFilePath);

  fp_r = fopen(outputFilePath, "w");
  if (datatype == DOUBLE)
    icount = fwrite(d_r, outSize, 1, fp_r);
  else /* float */
    icount = fwrite(f_r, outSize, 1, fp_r);
  if (icount != 1) {
    printf("Write the reconstructed file failed: %lu != %d!\n", outSize, icount);
    exit(1);
  }
  fclose(fp_r);

  int outliers=0;
  if (!strcmp(argv[1], "-d")) { /* double */
    double max_diff = 0.0;
    for (int i=0; i<N; i++) {
      double diff = fabs(d[i]-d_r[i]);
      if (diff > DBL_EPSILON) {
	if (outliers < 5)
          printf("reconstruction error=%e\n", diff);
        outliers++;
      }
      if (diff > max_diff)
	max_diff = diff;
    }
    printf("max_diff=%e\n", max_diff);
  }
  else {
    float max_diff = 0.0;
    for (int i=0; i<N; i++) {
      float diff = fabsf(f[i]-f_r[i]);
      if (diff > FLT_EPSILON) {
	if (outliers < 5)
          printf("reconstruction error=%e\n", diff);
        outliers++;
      }
      if (diff > max_diff)
	max_diff = diff;
    }
    printf("max_diff=%e\n", max_diff);
  }

  if (outliers != 0)
    printf("reconstructed data (in double) differ from the original (%d out of %d)\n", outliers, N);
    
  if (datatype == DOUBLE) { /* double */
    free(d); free(d_x); free(d_r);
  }
  else { /* float */
    free(f); free(f_x); free(f_r);
  }
  printf("done\n");

  return 0;
}
