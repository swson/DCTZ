/**
 * @file dctz-test.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief DCTZ test program for Z-Checker
 * (C) 2019 University of Massachusetts Lowell.
       See LICENSE in the top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dctz.h"
#ifdef WITH_Z_CHECKER
#include "zc.h"
#endif

int main(int argc, char * argv[])
{
  size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0;
  size_t type_size = 0;
  char *oriFilePath, outputFilePath[640];
#ifdef WITH_Z_CHECKER
  char *solName = NULL;
#endif
  char *varName;
  double error_bound;
  t_var *var_r; /* buffer for reconstructed data */
  t_var *var; /* buffer for original data */
  t_var *var_z; /* buffer for compressed data */
  int N, min_argc;
  t_datatype datatype;

#ifdef WITH_Z_CHECKER
  min_argc = 7;
#else
  min_argc = 6;
#endif

  if (argc < min_argc) {
#ifdef WITH_Z_CHECKER
    printf("Test case: %s -d|-f [err bound] [var name] [srcFilePath] [dimension sizes...] solName \n", argv[0]);
    printf("Example: %s -d 1E-3 sedov testdata/x86/testfloat_8_8_128.dat 8 8 128 dctz-ec(1E-3) \n", argv[0]);
#else
    printf("Test case: %s -d|-f [err bound] [var name] [srcFilePath] [dimension sizes...] \n", argv[0]);
    printf("Example: %s -d 1E-3 sedov testdata/x86/testfloat_8_8_128.dat 8 8 128 \n", argv[0]);
#endif
    exit(0);
  }

  error_bound = atof(argv[2]);
  varName = argv[3];

  assert (argc >= 6);

#ifdef WITH_Z_CHECKER
  if (argc >= 7) { /* 1D */
    r1 = N = atoi(argv[5]);
    solName = argv[6]; /* dummy when z-checker is not set */
  }
  if (argc >= 8) { /* 2D */
    r2 = atoi(argv[6]);
    N = r1*r2;
    solName = argv[7]; /* dummy when z-checker is not set */
  }
  if (argc >= 9) { /* 3D */
    r3 = atoi(argv[7]);
    N = r1*r2*r3;
    solName = argv[8]; /* dummy when z-checker is not set */
  }
  if (argc >= 10) { /* 4D */
    r4 = atoi(argv[8]);
    N = r1*r2*r3*r4;
    solName = argv[9]; /* dummy when z-checker is not set */
  }
#else
    if (argc >= 6) { /* 1D */
    r1 = N = atoi (argv[5]);
  }
  if (argc >= 7) { /* 2D */
    r2 = atoi(argv[6]);
    N = r1*r2;
  }
  if (argc >= 8) { /* 3D */
    r3 = atoi(argv[7]);
    N = r1*r2*r3;
  }
  if (argc >= 9) { /* 4D */
    r4 = atoi(argv[8]);
    N = r1*r2*r3*r4;
  }
#endif

  printf("total number of elements = %d\n", N);

  oriFilePath = argv[4];

  /* argv[2] = error bound, e.g., 1E-3 */
#ifdef USE_QTABLE
  sprintf(outputFilePath, "%s.qt.%s.z", oriFilePath, argv[2]);
#else
  sprintf(outputFilePath, "%s.ec.%s.z", oriFilePath, argv[2]);
#endif /* USE_QTABLE */

#ifdef WITH_Z_CHECKER
  ZC_Init("zc.config"); /* hard coded */
#endif /* WITH_Z_CHECKER */

  size_t outSize;
#ifdef WITH_Z_CHECKER
  ZC_DataProperty* dataProperty = NULL;
  ZC_CompareData *compareResult = NULL;
#endif /* WITH_Z_CHECKER */
  FILE *fp_in = fopen(oriFilePath, "rb");
  if (fp_in == NULL) {
    perror("Failed: ");
    printf("File Not Found\n");
    return (1);
  }

  if (!strcmp(argv[1], "-d")) {
    type_size = sizeof(double);
    datatype = DOUBLE;
  }
  else {
    type_size = sizeof(float);
    datatype = FLOAT;;
  }

  var = malloc(sizeof(t_var));
  var_r = malloc(sizeof(t_var));
  var_z = malloc(sizeof(t_var));

  if (datatype == DOUBLE) {
    if (NULL == (var->buf.d = (double *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: org_buf\n");
      exit(1);
    }
    if (NULL == (var_r->buf.d = (double *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: reconst_buf\n");
      exit(1);
    }
    if (NULL == (var_z->buf.d = (double *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: comp_buf\n");
      exit(1);
    }
    var->datatype = var_r->datatype = var_z->datatype = DOUBLE;
  }
  else { /*FLOAT */
    if (NULL == (var->buf.f = (float *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: org_buf\n");
      exit(1);
    }
    if (NULL == (var_r->buf.f = (float *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: reconst_buf\n");
      exit(1);
    }
    if (NULL == (var_z->buf.f = (float *)malloc(N*type_size))) {
      fprintf(stderr, "Out of memory: comp_buf\n");
      exit(1);
    }
    var->datatype = var_r->datatype = var_z->datatype = FLOAT;
  }

  size_t bytes_read;
  if (datatype == DOUBLE)
    bytes_read = fread(var->buf.d, type_size, N, fp_in);
  else /* FLOAT */
    bytes_read = fread(var->buf.f, type_size, N, fp_in);

  if (bytes_read != N) {
    perror("Error reading file");
    exit(EXIT_FAILURE);
  }
#ifdef WITH_Z_CHECKER
  if (datatype == DOUBLE)
    dataProperty = ZC_startCmpr(varName, ZC_DOUBLE, var->buf.d, r5, r4, r3, r2, r1);
  else /* FLOAT */
    dataProperty = ZC_startCmpr(varName, ZC_FLOAT, var->buf.f, r5, r4, r3, r2, r1);
#endif /* WITH_Z_CHECKER */
  dctz_compress(var, N, &outSize, var_z, error_bound);

  printf("oriFilePath = %s, outputFilePath = %s, datatype = %s error = %s, dim1 = %zu dim2 = %zu dim3 = %zu dim4 = %zu\n", oriFilePath, outputFilePath, datatype==FLOAT?"float":"double", argv[2], r1, r2, r3, r4);
  printf("outsize = %zu\n", outSize);

  /* restore original data for calculating PSNR properly */
  /* this is required for calc_psnr() and Z-checker */
  struct header h;
  if (datatype == DOUBLE) {
    memcpy(&h, var_z->buf.d, sizeof(struct header));
    /* deapply scaling factor */
    if (h.scaling_factor.d != 1.0) {
      int i;
      for (i=0; i<N; i++)
	var->buf.d[i] *= h.scaling_factor.d;
    }
  }
  else { /* FLOAT */
    memcpy(&h, var_z->buf.f, sizeof(struct header));
    /* deapply scaling factor */
    if (h.scaling_factor.f != 1.0) {
      int i;
      for (i=0; i<N; i++)
	var->buf.f[i] *= h.scaling_factor.f;
    }
  }

#ifdef WITH_Z_CHECKER
  compareResult = ZC_endCmpr(dataProperty, solName, outSize);
#endif /* WITH_Z_CHECKER */

  fclose(fp_in);

  char zfile[640];
  FILE *fp_z;
  int icount;

#ifdef USE_QTABLE
  sprintf(zfile, "%s.qt.%s.z", oriFilePath, argv[2]);
#else
  sprintf(zfile, "%s.ec.%s.z", oriFilePath, argv[2]);
#endif
  fp_z = fopen(zfile, "wb");

  if (datatype == DOUBLE)
    icount = fwrite(var_z->buf.d, outSize, 1, fp_z);
  else /* FLOAT */
    icount = fwrite(var_z->buf.f, outSize, 1, fp_z);
  if (icount != 1) {
    printf("Write qtz file failed: %lu != %d!\n", outSize, icount);
    exit(1);
  }
  fclose(fp_z);

#ifdef USE_QTABLE
  sprintf(zfile, "%s.qt.%s.z.r", oriFilePath, argv[2]);
#else
  sprintf(zfile, "%s.ec.%s.z.r", oriFilePath, argv[2]);
#endif /* USE_QTABLE */
  FILE *fp_r;
  fp_r = fopen(zfile, "wb");
#ifdef WITH_Z_CHECKER
  ZC_startDec();
#endif /* WITH_Z_CHECKER */

  dctz_decompress(var_z, var_r);
#ifdef WITH_Z_CHECKER
  if (datatype == DOUBLE)
    ZC_endDec(compareResult, var_r->buf.d);
  else /* FLOAT */
    ZC_endDec(compareResult, var_r->buf.f);
#endif /* WITH_Z_CHECKER */

  if (datatype == DOUBLE)
    icount = fwrite(var_r->buf.d, N*type_size, 1, fp_r);
  else /* FLOAT */
    icount = fwrite(var_r->buf.f, N*type_size, 1, fp_r);
  if (icount != 1) {
    printf("Write qtz.r file failed:  != %d!\n",  icount);
    exit(1);
  }

  fclose(fp_r);

#ifdef WITH_Z_CHECKER
  freeDataProperty(dataProperty);
  freeCompareResult(compareResult);
#endif /* WITH_Z_CHECKER */

  double cr, psnr;
  cr = (double)(N*type_size)/(double)outSize;
  psnr = calc_psnr(var, var_r, N, error_bound);
  printf("CR = %f, PSNR = %f\n", cr, psnr);

  free(var_z);
  free(var_r);
  free(var);
  printf("done\n");

#ifdef WITH_Z_CHECKER
  ZC_Finalize();
#endif /* WITH_Z_CHECKER */

  return 0;
}
