/**
 * @file dctz-dump.c
 * @author Seung Woo Son
 * @date March 2022
 * @brief DCTZ header dump program
 * (C) 2019-2022 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../dctz.h"

/* gcc -o dctz-dump dctz-dump.c -Wall -g */

int main(int argc, char * argv[])
{
  struct header h;
  unsigned char *buf;
  FILE *fp;

  if (argc != 2) {
    printf("Usage: %s filename\n", argv[0]);
    exit(0);
  }

  fp = fopen(argv[1], "rb");
  if (fp == NULL) {
    perror("Failed: ");
    printf("File Not Found\n");
    return 0;
  }

  buf = (unsigned char *)malloc(sizeof(struct header));

  size_t bytes_read;
  bytes_read = fread(buf, sizeof(h), 1, fp);

  memcpy(&h, buf, sizeof(struct header));

  printf("File Name=%s\n", argv[1]);
  printf("N=%d\n", h.num_elements);
  printf("error_bound=%f\n", h.error_bound);
  printf("total # of AC_exact=%d\n", h.tot_AC_exact_count);
  printf("SF=%f\n", h.scaling_factor.d);
  
  free(buf);
  fclose(fp);
  
  return 0;
}
