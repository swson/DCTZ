/**
 * @file dct.h
 * @author Seung Woo Son
 * @date July 2019
 * @brief header for DCT routines
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#ifndef _DCT_H_
#define _DCT_H_

#ifndef M_PI
# define M_PI 3.14159265358979323846 /* pi */
#endif

void dct_fftw(double *a, double *b, int dn, int nblk);
void dct_fftw_f(float *a, float *b, int dn, int nblk);
void ifft_idct(int dn, double *a, double *data);
void ifft_idct_f(int dn, float *a, float *data);
//void idct_fftw(void *a, void *b, int dn, int blk_i, int nblk, unsigned int datatype);
void dct_init(int dn);
void dct_init_f(int dn);
void dct_finish();
void dct_finish_f();
void idct_finish();
void idct_finish_f();

#endif
