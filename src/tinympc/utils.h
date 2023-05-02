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
#define PrintMatrix(mat)      \
  {                          \
    printf("%s = \n", #mat); \
    slap_PrintMatrix(mat);   \
  }

#define PrintMatrixT(mat)                   \
  {                                        \
    printf("%s = \n", #mat);               \
    slap_PrintMatrix(slap_Transpose(mat)); \
  }

//========================================
// Return length of an array
//========================================
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

//========================================
// Initialize memory with zeros
//========================================
#define INIT_ZEROS(data) (memset(data, 0, sizeof(data)))  

//========================================
// Return a random noise from percentage
//========================================
#define NOISE(percent) (((2 * ((float)rand() / RAND_MAX)) - 1) / 100 * percent)

//========================================
// Print matrix info
//========================================
#define PrintMatrixInfo(mat)      \
  {                          \
    printf("%s info: \n", #mat); \
    printf(" Dims: (%d, %d)\n", mat.rows, mat.cols);   \
    printf(" Data: "); \
    for (int imat = 0; imat < mat.cols * mat.rows; ++imat) { \
      printf("%.4f, ", mat.data[imat]); \
    } \
    printf("\n");\
  }

//========================================
// Read data from file
//========================================
int tiny_ReadData(const char* filename, sfloat* des, const int size,
                  bool verbose);

//========================================
// Read data from file and copy the last knot point into
// remaining space of the array. Useful for extend horizon at the end.
//========================================
int tiny_ReadData_Extend(const char* filename, sfloat* des, const int stride,
                         const int size, bool verbose);

//========================================
// Read data from file and copy the goal state into
// remaining space of the array. Useful for extend horizon at the end.
//========================================
int tiny_ReadData_ExtendGoal(const char* filename, sfloat* des,
                             const sfloat* xf, const int stride, const int size,
                             bool verbose);

//========================================
// Clamp the inputs to within min max value,
// will modify the provided array
//========================================
void tiny_Clamps(sfloat* arr, const sfloat* min, const sfloat* max,
                 const int N);

void tiny_Clamp(sfloat* arr, const sfloat min, const sfloat max, const int N);

void tiny_ClampMatrix(Matrix* mat, const Matrix min, const Matrix max);

void tiny_ShiftFill(Matrix* mats, const int length);

void tiny_ShiftFillWith(Matrix* mats, const sfloat* x, const int length);