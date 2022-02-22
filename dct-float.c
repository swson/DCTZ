/**
 * @file dct-float.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief DCT implementation for single precision data
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "dct.h"
#include "dctz.h"

static int flag = 0;
static fftwf_plan p = NULL;
static fftwf_complex *in = NULL, *out = NULL;
static float *as = NULL, *ax = NULL;
static float *ias = NULL, *iax = NULL;

void dct_init_f(int dn)
{
  int i;

  in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * 2*dn);
  out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * 2*dn);
  size_t typesize = sizeof(float);
  as = (float *)malloc(dn*typesize);
  ax = (float *)malloc(dn*typesize);

  float x = 0.0;
  float y;
  
  /* Compute weights to multiply DFT coefficients */
  for (i=0; i<dn; i++) {
    y = -i*M_PI/(2*dn);
    as[i] = exp(x)*cos(y)/sqrt(2.0*dn);
    ax[i] = exp(x)*sin(y)/sqrt(2.0*dn);
  }
  as[0] = as[0]/sqrt(2.0);
  if (!(dn%2)) {
    for (i=0; i<dn; i++) {
      as[i] = as[i]*2;
      ax[i] = ax[i]*2;
    }
    p = fftwf_plan_dft_1d(dn, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  }
  else {
    p = fftwf_plan_dft_1d(2*dn, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  }  
}

void dct_fftw_f(float *a, float *b, int dn,  int nblk) 
{
  int i, j, k;

  if (dn%2 == 1) { /* for odds */
    memset(in, 0, sizeof(fftwf_complex)*2*dn);
    for (i=0; i<dn; i++) {
      in[i][0] = a[i];
      in[dn+i][0] = a[dn-1-i];
    }
#ifdef DCT_DEBUG
    printf("\n the in out before execute: \n");
    for (i=0; i<dn; i++) {
      printf("in :%f  %f      \n ", in[i][0], in[i][1]);
      printf("out :%f  %f     \n ", out[i][0], out[i][1]);
    }
#endif
    fftwf_execute(p); /* repeat as needed*/
  } else { /*  for even */
    memset(in, 0, sizeof(fftwf_complex)*2*dn);
    for (i=0, j=0, k=dn-1; i<dn; i++) {
      if (i%2) {
	in[k][0] = a[i];
	--k;
      } else {
	in[j][0] = a[i];
	++j;
      }
    }
#ifdef DCT_DEBUG
    printf("\n the in out before execute: \n");
    for (i=0; i<dn; i++) {
      printf("in :%f  %f      \n ", in[i][0], in[i][1]);
      printf("out :%f  %f     \n ", out[i][0], out[i][1]);
    }
#endif
    fftwf_execute(p); /* repeat as needed*/
  }
#ifdef DCT_DEBUG
  printf("the in  after execute: \n");
  for (i=0; i<dn; i++) {
    printf("in :%f  %f      \n ", in[i][0], in[i][1]);
    printf("out :%f  %f     \n ", out[i][0], out[i][1]);
  }
#endif
  for (i=0; i<dn; i++) {
    b[i] = as[i]*out[i][0] - ax[i]*out[i][1];
  }
}

void dct_finish_f() {
  flag = 0;
  fftwf_destroy_plan(p);
  fftwf_free(in);
  fftwf_free(out);
  free(as);
  free(ax);
  fftwf_cleanup();
}

void ifft_idct_f(int dn, float *a,  float *data)
{
  int i, j, k;
  float ias_0;

  size_t typesize = sizeof(float);

  if (flag == 0) {
    in = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 2*dn); /* IFFT input */
    out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 2*dn); /* IFFT output */
    ias = (float *)malloc(dn*typesize);
    iax = (float *)malloc(dn*typesize);

    float x = 0.0;
    float y;

    /* Compute weights to multiply IDFT coefficients */
    for (i = 0; i < dn; i++) {
      y = i*M_PI/(2*dn);
      ias[i] = exp(x)*cos(y)*sqrt(2.0*dn);
      iax[i] = exp(x)*sin(y)*sqrt(2.0*dn);
    }
#ifdef DCT_DEBUG
    for (i=0; i<dn; i++) {
      printf("%f, %f\n", ias[i], iax[i]);
    }
#endif
  }
  memset(in, 0, sizeof(fftwf_complex)*dn*2);
  memset(out, 0, sizeof(fftwf_complex)*dn*2);
  
  if (dn%2 == 1) { /* for odd */
    ias_0 = ias[0] * sqrt(2.0);
    in[0][0] = ias_0*a[0];
    in[0][1] = iax[0]*a[0];

    for(i = 1; i < dn; i++) {
      in[i][0] = ias[i]*a[i];
      in[i][1] = iax[i]*a[i];
      in[dn+i][0] = iax[i]*a[dn-i];
      in[dn+i][1] = -ias[i]*a[dn-i];
    }

    if (flag == 0) {
      p = fftwf_plan_dft_1d(2*dn, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
      flag = 1;
    }
    fftwf_execute(p); /* repeat as needed*/
    
    for (i=0; i<dn; i++) {
      data[i] = out[i][0] / dn / 2;
    }
  } else { /* for even */
    ias_0 = ias[0] / sqrt(2.0);
    in[0][0] = ias_0*a[0];
    in[0][1] = iax[0]*a[0];

    for (i=1; i<dn; i++) {
      in[i][0] = ias[i]*a[i];
      in[i][1] = iax[i]*a[i];
#ifdef DCT_DEBUG
      printf("%f,  %f\n", in[i][0], in[i][1]);
#endif
    }
    
    if (flag == 0) {
      p = fftwf_plan_dft_1d(dn, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
      flag = 1;
    }
    fftwf_execute(p); /* repeat as needed*/
    
    for (i=0; i<dn; i++) {
      out[i][0] = out[i][0] / dn;
      out[i][1] = out[i][1] / dn;
#ifdef DCT_DEBUG
      printf("out == %f\n", out[i][0]);
#endif
    }
    
    for (i=0, j=0, k=dn-1; i<dn; i++) {
      if (i%2) {
	data[i] = out[k][0];
	k--;
      } else{
	data[i] = out[j][0];
	++j;
      }
#ifdef DCT_DEBUG
      printf("data ==  %f\n", data[i]);
#endif
    }
  }
}

void idct_finish_f()
{
  flag = 0;
  fftwf_destroy_plan(p);
  fftwf_free(in);
  fftwf_free(out);
  free(iax);
  free(ias);
  fftwf_cleanup();
}
