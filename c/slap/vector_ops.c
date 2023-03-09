//
// Created by Brian Jackson on 12/18/22.
// Copyright (c) 2022 Robotic Exploration Lab. All rights reserved.
//

#include "vector_ops.h"

#include <math.h>

#include "matrix_checks.h"

MatrixIterator slap_ArgMax(Matrix mat, double* max_value) {
  MatrixIterator max_index = slap_Iterator(mat);
  double value = -INFINITY;
  double value_i;
  for (MatrixIterator it = slap_Iterator(mat); !slap_IsFinished(&it);
       slap_Step(&it)) {
    value_i = mat.data[it.index];
    if (value_i > value) {
      value = value_i;
      max_index = it;
    }
  }
  if (max_value) *max_value = value;
  return max_index;
}

double slap_Max(Matrix mat) {
  double max_value;
  slap_ArgMax(mat, &max_value);
  return max_value;
}

MatrixIterator slap_ArgMin(Matrix mat, double* min_value) {
  MatrixIterator min_index = slap_Iterator(mat);
  double value = +INFINITY;
  double value_i;
  for (MatrixIterator it = slap_Iterator(mat); !slap_IsFinished(&it);
       slap_Step(&it)) {
    value_i = mat.data[it.index];
    if (value_i < value) {
      value = value_i;
      min_index = it;
    }
  }
  if (min_value) *min_value = value;
  return min_index;
}

double slap_Min(Matrix mat) {
  double min_value;
  slap_ArgMin(mat, &min_value);
  return min_value;
}

double slap_NormTwoSquared(Matrix mat) {
  SLAP_ASSERT_DENSE(mat, NAN,
                    "NormTwoSquared must be called on a dense matrix");
  double value = 0;
  double value_i;
  for (int i = 0; i < slap_NumElements(mat); ++i) {
    value_i = mat.data[i];
    value += value_i * value_i;
  }
  return value;
}

double slap_NormTwo(Matrix mat) {
  SLAP_ASSERT_DENSE(mat, NAN, "NormTwo must be called on a dense matrix");
  double norm_squared = slap_NormTwoSquared(mat);
  return sqrt(norm_squared);
}

double slap_NormInf(Matrix mat) {
  SLAP_ASSERT_DENSE(mat, NAN, "NormInf must be called on a dense matrix");
  double value = 0;
  double value_i;
  for (int i = 0; i < slap_NumElements(mat); ++i) {
    value_i = fabs(mat.data[i]);
    if (value_i > value) {
      value = value_i;
    }
  }
  return value;
}

double slap_NormOne(Matrix mat) {
  SLAP_ASSERT_DENSE(mat, NAN, "NormOne must be called on a dense matrix");
  double value = 0;
  double value_i;
  for (int i = 0; i < slap_NumElements(mat); ++i) {
    value_i = fabs(mat.data[i]);
    value += value_i;
  }
  return value;
}

double slap_Sum(Matrix mat) {
  double sum = 0;
  for (MatrixIterator it = slap_Iterator(mat); !slap_IsFinished(&it);
       slap_Step(&it)) {
    double value_i = mat.data[it.index];
    sum += value_i;
  }
  return sum;
}
