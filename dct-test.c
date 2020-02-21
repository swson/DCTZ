/**
 * @file dct-test.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief testing DCT performance
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dctz.h"

/* gcc -o dct-test dct-test.c dct.c dct-float.c -lfftw3 -lfftw3-f -lm -Wall -g */

enum numtype {DOUBLE, FLOAT};

int main (int argc, char * argv[])
{
  size_t r4=0,r3=0,r2=0,r1=0;
  size_t typesize = 0;
  char oriFilePath[640], outputFilePath[640];
  double *d, *d_x;
  float *f, *f_x;
  int datatype;
  int N, i, nblk;;
  
  if (argc < 4) {
    printf ("Test case: %s -d|-f [srcFilePath] [dimension sizes...]\n", argv[0]);
    printf ("Example: %s -d testdata/x86/testfloat_8_8_128.dat 8 8 128 dctz-ec(1E-3) \n", argv[0]);

    exit (0);
  }
  
  assert (argc >= 4);

  if (argc >= 4) { /* 1D */
    r1 = N = atoi (argv[3]);
  }
  if (argc >= 5) { /* 2D */
    r2 = atoi (argv[4]);
    N = r1*r2;
  }
  if (argc >= 6) { /* 3D */
    r3 = atoi (argv[5]);
    N = r1*r2*r3;
  } 
  if (argc >= 7) { /* 4D */
    r4 = atoi (argv[6]);
    N = r1*r2*r3*r4;
  }
  
  printf ("total number = %d\n", N);

  nblk = ceil(N/BLK_SZ);

  sprintf (oriFilePath, "%s", argv[2]);
  
  FILE *fp_in = fopen (oriFilePath, "rb");
  if (fp_in == NULL) {
    perror ("Failed: ");
    printf ("File Not Found\n");
    return (1);
  }
  
  if (!strcmp (argv[1], "-d")) { /* double */
    typesize = sizeof(double);
    datatype = DOUBLE;
    if (NULL == (d = (double *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a\n");
      exit (1);
    }

    if (NULL == (d_x = (double *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a_z\n");
      exit (1);
    }
    fread (d, typesize, N, fp_in);
    dct_init (BLK_SZ);
    for (i=0; i<nblk; i++)
      dct_fftw (d+i*BLK_SZ, d_x+i*BLK_SZ, BLK_SZ, nblk); // use min(BLK_SZ,)
  } 
  else { /* float */
    typesize = sizeof (float);
    datatype = FLOAT;
    if (NULL == (f = (float *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a\n");
      exit (1);
    }
    if (NULL == (f_x = (float *)malloc (N*typesize))) {
      fprintf(stderr, "Out of memory: a\n");
      exit (1);
    }
    fread (f, typesize, N, fp_in);
    dct_init_f (BLK_SZ);
    for (i=0; i<nblk; i++)
      dct_fftw_f (f+i*BLK_SZ, f_x+i*BLK_SZ, BLK_SZ, nblk); // use min(BLK_SZ,)
  }	  
  
  fclose (fp_in);

  FILE *fp_x;
  int icount;
  size_t outSize = N; /* outSize is same as N */
  
  sprintf (outputFilePath, "%s.x", oriFilePath);

  fp_x = fopen (outputFilePath, "wb");
  if (datatype == DOUBLE)
    icount = fwrite (d_x, outSize, 1, fp_x);
  else /* float */
    icount = fwrite (f_x, outSize, 1, fp_x);
  if (icount != 1) {
    printf ("Write DCT coefficients to a file failed: %lu != %d!\n", outSize, icount);
    exit (1);
  }
  fclose (fp_x);

  if (datatype == DOUBLE) { /* double */
    free (d); free (d_x);
  }
  else { /* float */
    free (f); free (f_x);
  }
  printf ("done\n");

  return 0;
}
