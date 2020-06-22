/**
 * @file dctz-test.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief DCTZ test program for Z-Checker
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dctz.h"
#ifdef WITH_Z_CHECKER
#include "zc.h"
#endif

int main (int argc, char * argv[])
{
  size_t r5=0,r4=0,r3=0,r2=0,r1=0;
  size_t typesize = 0;
  char oriFilePath[640], outputFilePath[640];
#ifdef WITH_Z_CHECKER
  char *solName = NULL;
#endif
  char *varName; 
  double error_bound;
  void *a_r;
  double *d;
  float *f;
  int datatype;
  char *a_z;
  int N, min_argc;

#ifdef WITH_Z_CHECKER
  min_argc = 7;
#else
  min_argc = 6;
#endif
  
  if (argc < min_argc) {
#ifdef WITH_Z_CHECKR
    printf ("Test case: %s -d|-f [err bound] [var name] [srcFilePath] [dimension sizes...] solName \n", argv[0]);
    printf ("Example: %s -d 1E-3 sedov testdata/x86/testfloat_8_8_128.dat 8 8 128 dctz-ec(1E-3) \n", argv[0]);
#else
    printf ("Test case: %s -d|-f [err bound] [var name] [srcFilePath] [dimension sizes...] \n", argv[0]);
    printf ("Example: %s -d 1E-3 sedov testdata/x86/testfloat_8_8_128.dat 8 8 128 \n", argv[0]);
#endif
    exit (0);
  }
  
  error_bound = atof (argv[2]);
  varName=argv[3];

  assert (argc >= 6);

#ifdef WITH_Z_CHECKER
  if (argc >= 7) { /* 1D */
    r1 = N = atoi (argv[5]);
    solName = argv[6]; /* dummy when z-checker is not set */
  }
  if (argc >= 8) { /* 2D */
    r2 = atoi (argv[6]);
    N = r1*r2;
    solName = argv[7]; /* dummy when z-checker is not set */
  }
  if (argc >= 9) { /* 3D */
    r3 = atoi (argv[7]);
    N = r1*r2*r3;
    solName = argv[8]; /* dummy when z-checker is not set */
  } 
  if (argc >= 10) { /* 4D */
    r4 = atoi (argv[8]);
    N = r1*r2*r3*r4;
    solName = argv[9]; /* dummy when z-checker is not set */
  }
#else
    if (argc >= 6) { /* 1D */
    r1 = N = atoi (argv[5]);
  }
  if (argc >= 7) { /* 2D */
    r2 = atoi (argv[6]);
    N = r1*r2;
  }
  if (argc >= 8) { /* 3D */
    r3 = atoi (argv[7]);
    N = r1*r2*r3;
  } 
  if (argc >= 9) { /* 4D */
    r4 = atoi (argv[8]);
    N = r1*r2*r3*r4;
  }
#endif
  
  printf ("total number = %d\n", N);
	
  sprintf (oriFilePath, "%s", argv[4]);
  
#ifdef USE_QTABLE
  sprintf (outputFilePath, "%s.qt.%s.z", oriFilePath, argv[2]);
#else
  sprintf (outputFilePath, "%s.t.%s.z", oriFilePath, argv[2]);
#endif /* USE_QTABLE */

#ifdef WITH_Z_CHECKER
  ZC_Init ("zc.config"); /* hard coded */
#endif /* WITH_Z_CHECKER */
  
  size_t outSize;
#ifdef WITH_Z_CHECKER
  ZC_DataProperty* dataProperty = NULL;
  ZC_CompareData *compareResult = NULL;
#endif /* WITH_Z_CHECKER */
  FILE *fp_in = fopen (oriFilePath, "rb");
  if (fp_in == NULL) {
    perror ("Failed: ");
    printf ("File Not Found\n");
    return (1);
  }
  
  if (!strcmp (argv[1], "-d")) {
    typesize = sizeof(double);
    datatype = data_type_double;
    if (NULL == (d = (double *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a\n");
      exit (1);
    }
    if (NULL == (a_r = (double *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a\n");
      exit (1);
    }
    if (NULL == (a_z = (char *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a_z\n");
      exit (1);
    }

    size_t bytes_read = fread (d, typesize, N, fp_in);
    if (bytes_read != N) {
      perror ("Error reading file");
      exit (EXIT_FAILURE);
    }
#ifdef WITH_Z_CHECKER
    dataProperty = ZC_startCmpr (varName, ZC_DOUBLE, d, r5, r4, r3, r2, r1);
#endif /* WITH_Z_CHECKER */
    dctz_compress (d, N, &outSize, a_z, error_bound);
  } 
  else {	
    typesize = sizeof (float);
    datatype = data_type_float;
    if (NULL == (f = (float *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a\n");
      exit (1);
    }
    if (NULL == (a_r = (float *)malloc (N*typesize))) {
      fprintf(stderr, "Out of memory: a\n");
      exit (1);
    }
    if (NULL == (a_z = (char *)malloc (N*typesize))) {
      fprintf (stderr, "Out of memory: a_z\n");
      exit (1);
    }
    size_t bytes_read = fread (f, typesize, N, fp_in);
    if (bytes_read != N) {
      perror ("Error reading file");
      exit (EXIT_FAILURE);
    }
#ifdef WITH_Z_CHECKER
    dataProperty = ZC_startCmpr (varName, ZC_FLOAT, f, r5, r4, r3, r2, r1);
#endif /* WITH_Z_CHECKER */
    dctz_compress_float (f, N, &outSize, a_z, error_bound);
  }	  
  printf ("oriFilePath = %s, outputFilePath = %s, datatype = %s error = %s, dim1 = %zu dim2 = %zu dim3 = %zu dim4 = %zu\n", oriFilePath, outputFilePath, datatype==0?"double":"float", argv[2], r1, r2, r3, r4);
  printf ("outsize = %zu\n", outSize);
#ifdef WITH_Z_CHECKER
  compareResult = ZC_endCmpr (dataProperty, solName, outSize);
#endif /* WITH_Z_CHECKER */
  
  struct header h;
  
  memcpy (&h, a_z, sizeof(struct header));
  double SF = h.scaling_factor;
  
#ifdef DEBUG
  printf ("SF = %f\n", SF);
#endif /* DEBUG */ 
  // deapply scaling factor to the original data
  double xscale = pow (10, SF-1);
  if (SF != 1.0)
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(a, SF)
#endif
    for (int i=0; i<N; i++) {
      if (datatype == data_type_double) 
	d[i] *= xscale; 
      else 
	f[i] *= xscale;
    }
#ifdef DEBUG
  for (int i=0; i<BLK_SZ; i++) { // show the first block
    printf ("d[%d] = %e %p\n", i, d[i], &d[i]);
    if (i%BLK_SZ == 0 && i != 0) printf ("\n");
  }
#endif
  fclose (fp_in);
  
  char zfile[640];
  FILE *fp_z;
  int icount;
  
#ifdef USE_QTABLE
  sprintf (zfile, "%s.qt.%s.z", oriFilePath, argv[2]);
#else
  sprintf (zfile, "%s.t.%s.z", oriFilePath, argv[2]);
#endif
  fp_z = fopen (zfile, "wb");
  icount = fwrite (a_z, outSize, 1, fp_z);
  if (icount != 1) {
    printf ("Write qtz file failed: %lu != %d!\n", outSize, icount);
    exit (1);
  }
  fclose (fp_z);
  
#ifdef USE_QTABLE
  sprintf (zfile, "%s.qt.%s.z.r", oriFilePath, argv[2]);
#else
  sprintf (zfile, "%s.t.%s.z.r", oriFilePath, argv[2]);
#endif /* USE_QTABLE */
  FILE *fp_r;
  fp_r = fopen (zfile, "wb");
#ifdef WITH_Z_CHECKER
  ZC_startDec ();
#endif /* WITH_Z_CHECKER */
  
  if (datatype == data_type_double) {
    dctz_decompress (a_z, (double *) a_r);
#ifdef WITH_Z_CHECKER
    ZC_endDec (compareResult, (double *) a_r);
#endif /* WITH_Z_CHECKER */
    icount = fwrite ((double *)a_r, N*sizeof(double), 1, fp_r);
  }
  else {
    dctz_decompress_float (a_z, (float *) a_r);
#ifdef WITH_Z_CHECKER
    ZC_endDec (compareResult, (float *)a_r);
#endif /* WITH_Z_CHECKER */
    icount = fwrite ((float *)a_r, N*sizeof(float), 1, fp_r);
  }
  
  if (icount != 1) {
    printf ("Write qtz.r file failed:  != %d!\n",  icount);
    exit (1);
  }

  fclose (fp_r);

#ifdef WITH_Z_CHECKER
  freeDataProperty (dataProperty);
  freeCompareResult (compareResult);
#endif /* WITH_Z_CHECKER */
  free (a_z);
  free (a_r);
  printf ("done\n");

#ifdef WITH_Z_CHECKER
  ZC_Finalize ();
#endif /* WITH_Z_CHECKER */
  
  return 0;
}
