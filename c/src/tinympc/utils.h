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
#define tiny_Print(mat)      \
  {                          \
    printf("%s = \n", #mat); \
    slap_PrintMatrix(mat);   \
  }

#define tiny_PrintT(mat)                   \
  {                                        \
    printf("%s = \n", #mat);               \
    slap_PrintMatrix(slap_Transpose(mat)); \
  }

//========================================
// Return length of an array
//========================================
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

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