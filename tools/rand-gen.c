/**
 * @file rand-gen.c
 * @author Seung Woo Son
 * @date May 2022
 * @brief a program to generate a file with random numbers in range.
 * (C) 2019-2022 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main()
{
  int lower = 0, upper = 127, count = 10000; /* {0,127}, {128,256}, or (0,255} */
  int *buffer;
  FILE *fp;
  char fname[100];

  srand(time(0));

  buffer = malloc(count*sizeof(int));

  printf("Printing the first 10 randomly generated numbers... ");
  for (int i = 0; i < count; i++) {
    int num = (rand() % (upper - lower + 1)) + lower;
    if (i < 10) printf("%d ", num);
    buffer[i] = num;
  }

  printf("\n");

  sprintf(fname, "rand_%d_%d_%d.bin", lower, upper, count);
  fp = fopen(fname, "w");
  fwrite(buffer, 1, count*sizeof(int), fp);
  fclose(fp);
  
  return 0;
}
