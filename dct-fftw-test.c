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

int main(int argc, char *argv[])
{
  size_t r3=0, r2=0, r1=0;
  char *inFilePath, outFilePath[640];
  FILE *fin, *fout;
  double *d, *dx, *dr;
  float *f, *fx, *fr;
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
  
  fin = fopen(inFilePath, "r");

  if (!strcmp(argv[1], "-d")) { /* double */
    d = malloc(sizeof(double)*N);
    dx = malloc(sizeof(double)*N);
    dr = malloc(sizeof(double)*N);
    fread(d, sizeof(double), N, fin);
    fftw_plan plan = fftw_plan_r2r_1d(N, d, dx, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_plan plan_inv = fftw_plan_r2r_1d(N, dx, dr, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(plan_inv);

    /* FFTW doesn't rescale the output of the transform,
     * so the result must be divided by 2*N
     */
    for (int i=0; i<N; i++)
      dr[i] /= 2.0*N;
      
    fout = fopen(outFilePath, "w");
    fwrite(dr, sizeof(double)*N, 1, fout);

    int outliers=0;
    for (int i=0; i<N; i++) {
      if ((d[i]-dr[i]) > DBL_EPSILON) {
	if (outliers < 5)
	  printf("reconstruction error=%e\n", d[i]-dr[i]);
	outliers++;
      }
    }

    if (outliers != 0)
      printf("reconstructed data (in double) differ from the original (%d out of %d)\n", outliers, N);
    
    free(d); free(dx); free(dr);
  }
  else { /* float */
    f = malloc(sizeof(float)*N);
    fx = malloc(sizeof(float)*N);
    fr = malloc(sizeof(float)*N);
    fread(f, sizeof(float), N, fin);
    fftwf_plan plan = fftwf_plan_r2r_1d(N, f, fx, FFTW_REDFT10, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_plan plan_inv = fftwf_plan_r2r_1d(N, fx, fr, FFTW_REDFT01, FFTW_ESTIMATE);
    fftwf_execute(plan_inv);

    for (int i=0; i<N; i++)
      fr[i] /= 2.0f*N;

    fout = fopen(outFilePath, "w");
    fwrite(fr, sizeof(float)*N, 1, fout);

    int outliers=0;
    for (int i=0; i<N; i++) {
      if ((f[i]-fr[i]) > FLT_EPSILON) {
	if (outliers < 5)
	  printf("reconstruction error=%e\n", f[i]-fr[i]);
	outliers++;
      }
    }

    if (outliers != 0)
      printf("reconstructed data (in float) differ from the original (%d out of %d)\n", outliers, N);
    
    free(f); free(fx); free(fr);
  }

  fclose(fin);
  fclose(fout);

  /* Finding the differences between two binary files
     $ cmp fileA fileB
     or
     $ diff fileA fileB
   */
  
  return 0;
}
