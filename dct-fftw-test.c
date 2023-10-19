/**
 * @file dct-fftw-test.c
 * @author Seung Woo Son
 * @date April 2023
 * @brief DCT routine using fftw3
 * (C) 2023 University of Massachusetts Lowell.
       See LICENSE in the top-level directory.
*/

// gcc -o dct-fftw-test dct-fftw-test.c -lm -lfftw3 -lfftw3f
// sample dataset:
//   double-precision: https://sites.uml.edu/seungwoo-son/files/2019/07/dctz-test-data.zip
//   single-precision: http://www.mcs.anl.gov/~shdi/download/CESM-ATM-tylor.tar.gz

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <fftw3.h>
#include <math.h>

int main(int argc, char *argv[])
{
  size_t r3=0, r2=0, r1=0;
  char *inFilePath, xFilePath[640], outFilePath[640];
  FILE *fileIn, *fileX, *fileOut;
  double *d, *d_x, *d_r;
  float *f, *f_x, *f_r;
  int N;

  if (argc < 4) {
    printf("Test case: %s -d|-f [srcFilePath] [dimension size]\n", argv[0]);
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

  inFilePath = argv[2];
  sprintf(outFilePath, "%s.r", inFilePath);
  sprintf(xFilePath, "%s.x", inFilePath);
  
  fileIn = fopen(inFilePath, "r");

  if (!strcmp(argv[1], "-d")) { /**************** double ****************/
    d = malloc(sizeof(double)*N);
    d_x = malloc(sizeof(double)*N);
    d_r = malloc(sizeof(double)*N);
    fread(d, sizeof(double), N, fileIn);
    if (argc == 4) { /* 1D */
      fftw_plan plan = fftw_plan_r2r_1d(N, d, d_x, FFTW_REDFT10, FFTW_ESTIMATE);
      fftw_execute(plan);
      fftw_plan plan_inv = fftw_plan_r2r_1d(N, d_x, d_r, FFTW_REDFT01, FFTW_ESTIMATE);
      fftw_execute(plan_inv);

      /* FFTW doesn't rescale the output of the transform,
       * so the result must be divided by 2*N
       */
      for (int i=0; i<N; i++)
	d_r[i] /= 2.0*N;
    }
    else if (argc == 5) { /* 2D */
      fftw_plan plan = fftw_plan_r2r_2d(r1, r2, d, d_x, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
      fftw_execute(plan);
      fftw_plan plan_inv = fftw_plan_r2r_2d(r1, r2, d_x, d_r, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
      fftw_execute(plan_inv);

      /* FFTW doesn't rescale the output of the transform,
       * so the result must be divided by 2*N
       */
      for (int i=0; i<N; i++)
	d_r[i] /= 4.0*N;

    }
    else if (argc == 6) { /* 3D */
      fftw_plan plan = fftw_plan_r2r_3d(r1, r2, r3, d, d_x, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
      fftw_execute(plan);
      fftw_plan plan_inv = fftw_plan_r2r_3d(r1, r2, r3, d_x, d_r, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
      fftw_execute(plan_inv);

      /* FFTW doesn't rescale the output of the transform,
       * so the result must be divided by 2*N
       */
      for (int i=0; i<N; i++)
	d_r[i] /= 8.0*N;
    }

    fileX = fopen(xFilePath, "w");
    fwrite(d_x, sizeof(double)*N, 1, fileX);
    
    fileOut = fopen(outFilePath, "w");
    fwrite(d_r, sizeof(double)*N, 1, fileOut);

    int outliers=0;
    double max_diff=0.0;
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
    
    if (outliers != 0)
      printf("reconstructed data (in double) differ from the original (%d out of %d)\n", outliers, N);
    
    free(d); free(d_x); free(d_r);
  }
  else { /**************** float ****************/
    f = malloc(sizeof(float)*N);
    f_x = malloc(sizeof(float)*N);
    f_r = malloc(sizeof(float)*N);
    fread(f, sizeof(float), N, fileIn);

    if (argc == 4) { /* 1D */
      fftwf_plan plan = fftwf_plan_r2r_1d(N, f, f_x, FFTW_REDFT10, FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_plan plan_inv = fftwf_plan_r2r_1d(N, f_x, f_r, FFTW_REDFT01, FFTW_ESTIMATE);
      fftwf_execute(plan_inv);

      /* FFTW doesn't rescale the output of the transform,
       * so the result must be divided by 2*N
       */
      for (int i=0; i<N; i++)
	f_r[i] /= 2.0f*N;
    }
    else if (argc == 6) { /* 2D */
      fftwf_plan plan = fftwf_plan_r2r_2d(r1, r2, f, f_x, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_plan plan_inv = fftwf_plan_r2r_2d(r1, r2, f_x, f_r, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
      fftwf_execute(plan_inv);

      /* FFTW doesn't rescale the output of the transform,
       * so the result must be divided by 2^2*N
       */
      for (int i=0; i<N; i++)
	f_r[i] /= 4.0f*N;
    }
    else if (argc == 6) { /* 3D */
      fftwf_plan plan = fftwf_plan_r2r_3d(r1, r2, r3, f, f_x, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_plan plan_inv = fftwf_plan_r2r_3d(r1, r2, r3, f_x, f_r, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
      fftwf_execute(plan_inv);

      /* FFTW doesn't rescale the output of the transform,
       * so the result must be divided by 2^3*N
       */
      for (int i=0; i<N; i++)
	f_r[i] /= 8.0f*N;
    }

    fileX = fopen(xFilePath, "w");
    fwrite(f_x, sizeof(float)*N, 1, fileX);
    
    fileOut = fopen(outFilePath, "w");
    fwrite(f_r, sizeof(float)*N, 1, fileOut);

    int outliers=0;
    float max_diff=0.0;
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

    if (outliers != 0)
      printf("reconstructed data (in float) differ from the original (%d out of %d)\n", outliers, N);
    
    free(f); free(f_x); free(f_r);
  }

  fclose(fileIn);
  fclose(fileX);
  fclose(fileOut);

  /* Finding the differences between two binary files
     $ cmp fileA fileB
     or
     $ diff fileA fileB
   */
  
  return 0;
}
