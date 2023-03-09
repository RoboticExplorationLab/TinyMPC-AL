#pragma once

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "slap/slap.h"

//========================================
// Print matrix with its name (dummy)
//========================================
#define tiny_Print(mat)        \
  {                            \
    printf("\n%s = \n", #mat); \
    slap_PrintMatrix(mat);     \
  }

//========================================
// Read data from file
//========================================
int tiny_ReadData(const char* filename, double* des, const int size,
                  bool verbose) {
  FILE* input;
  int i;

  input = fopen(filename, "r");
  if (!input) {
    if (verbose == true)
      fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
    return EXIT_FAILURE;
  }

  for (i = 0; i < size; ++i) {
    if (fscanf(input, "%lf ", &(des[i])) != 1) {
      if (verbose == true) fprintf(stderr, "Invalid data in %s.\n", filename);
      fclose(input);
      return EXIT_FAILURE;
    }

    if (verbose == true) printf("Read %lf from %s.\n", des[i], filename);
  }

  if (ferror(input)) {
    fclose(input);
    if (verbose == true) fprintf(stderr, "Error reading %s.\n", filename);
    return EXIT_FAILURE;
  }
  if (fclose(input)) {
    if (verbose == true) fprintf(stderr, "Error closing %s.\n", filename);
    return EXIT_FAILURE;
  }

  if (verbose == true) printf("All doubles read successfully.\n");

  return EXIT_SUCCESS;
}

//========================================
// Read data from file and copy the last knot point into
// remaining space of the array. Useful for extend horizon at the end.
//========================================
int tiny_ReadData_Extend(const char* filename, double* des, const int stride,
                         const int size, bool verbose) {
  FILE* input;
  int i;
  int k = 0;
  input = fopen(filename, "r");
  if (!input) {
    if (verbose == true)
      fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
    return EXIT_FAILURE;
  }

  for (i = 0; i < size; ++i) {
    if (fscanf(input, "%lf ", &(des[i])) != 1) {
      if (verbose == true) fprintf(stderr, "Invalid data in %s.\n", filename);
      fclose(input);
      break;
    }

    if (verbose == true) printf("Read %lf from %s.\n", des[i], filename);

    k += 1;
  }

  if (verbose == true)
    printf("All doubles read successfully and now extend.\n");

  int remain_cnt = (size - k) / stride;  // # of remaining chunks
  for (i = 0; i < remain_cnt; i += 1) {
    for (int j = 0; j < stride; j += 1) {
      des[k + j + i * stride] = des[k + j - stride];  // copy
    }
  }

  return EXIT_SUCCESS;
}

//========================================
// Read data from file and copy the goal state into
// remaining space of the array. Useful for extend horizon at the end.
//========================================
int tiny_ReadData_ExtendGoal(const char* filename, double* des,
                             const double* xf, const int stride, const int size,
                             bool verbose) {
  FILE* input;
  int i;
  int k = 0;
  input = fopen(filename, "r");
  if (!input) {
    if (verbose == true)
      fprintf(stderr, "Cannot open %s: %s.\n", filename, strerror(errno));
    return EXIT_FAILURE;
  }

  for (i = 0; i < size; ++i) {
    if (fscanf(input, "%lf ", &(des[i])) != 1) {
      if (verbose == true) fprintf(stderr, "Invalid data in %s.\n", filename);
      fclose(input);
      break;
    }

    if (verbose == true) printf("Read %lf from %s.\n", des[i], filename);

    k += 1;
  }

  if (verbose == true)
    printf("All doubles read successfully and now extend.\n");

  int remain_cnt = (size - k) / stride;  // # of remaining chunks
  for (i = 0; i < remain_cnt; i += 1) {
    for (int j = 0; j < stride; j += 1) {
      des[k + j + i * stride] = xf[j];  // copy
    }
  }

  return EXIT_SUCCESS;
}

//========================================
// Clamp the inputs to within min max value
//========================================
void tiny_Clamps(double* arr, const double* min, const double* max,
                 const int N) {
  for (int k = 0; k < N; ++k) {
    arr[k] = (arr[k] > max[k]) ? max[k] : ((arr[k] < min[k]) ? min[k] : arr[k]);
  }
}
void tiny_Clamp(double* arr, const double min, const double max, const int N) {
  for (int k = 0; k < N; ++k) {
    arr[k] = (arr[k] > max) ? max : ((arr[k] < min) ? min : arr[k]);
  }
}