#include "matrix.h"

#ifndef PRECISION
#define PRECISION 5
#endif

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Matrix slap_NewMatrix(int rows, int cols) {
  double* data = (double*)malloc(rows * cols * sizeof(double));
  Matrix mat = {rows, cols, data};
  return mat;
}

Matrix slap_MatrixFromArray(int rows, int cols, double* data) {
  Matrix mat = {rows, cols, data};
  return mat;
}

Matrix slap_NewMatrixZeros(int rows, int cols) {
  double* data = (double*)calloc(rows * cols, sizeof(double));
  Matrix mat = {rows, cols, data};
  return mat;
}

int slap_MatrixSetConst(Matrix* mat, double val) {
  if (!mat) {
    return -1;
  }
  for (int i = 0; i < slap_MatrixNumElements(mat); ++i) {
    mat->data[i] = val;
  }
  return 0;
}

int slap_FreeMatrix(Matrix* mat) {
  if (mat) {
    if (mat->data) {
      free(mat->data);
      mat->data = NULL;
      return 0;
    }
  }
  return -1;
}

int slap_MatrixNumElements(const Matrix* mat) {
  if (!mat) {
    return -1;
  }
  return mat->rows * mat->cols;
}

int slap_MatrixGetLinearIndex(const Matrix* mat, int row, int col) {
  if (!mat) {
    return -1;
  }
  if (row < 0 || col < 0) {
    return -1;
  }
  return row + mat->rows * col;
}

double* slap_MatrixGetElement(Matrix* mat, int row, int col) {
  if (!mat) {
    return NULL;
  }
  return mat->data + slap_MatrixGetLinearIndex(mat, row, col);
}

const double* slap_MatrixGetElementConst(const Matrix* mat, int row, int col) {
  if (!mat) {
    return NULL;
  }
  return mat->data + slap_MatrixGetLinearIndex(mat, row, col);
}

double* slap_MatrixGetElementTranspose(Matrix* mat, int row, int col,
                                       bool istranposed) {
  double* out;
  if (!istranposed) {
    out = slap_MatrixGetElement(mat, row, col);
  } else {
    int row_transpose = col;
    int col_transpose = row;
    out = slap_MatrixGetElement(mat, row_transpose, col_transpose);
  }
  return out;
}

const double* slap_MatrixGetElementTransposeConst(const Matrix* mat, int row, int col,
                                       bool istranposed) {
  const double* out;
  if (!istranposed) {
    out = slap_MatrixGetElementConst(mat, row, col);
  } else {
    int row_transpose = col;
    int col_transpose = row;
    out = slap_MatrixGetElementConst(mat, row_transpose, col_transpose);
  }
  return out;
}

int slap_MatrixSetElement(Matrix* mat, int row, int col, double val) {
  if (!mat) {
    return -1;
  }
  int linear_index = slap_MatrixGetLinearIndex(mat, row, col);
  if (linear_index >= 0) {
    mat->data[linear_index] = val;
  } else {
    return -1;
  }
  return 0;
}

int slap_MatrixCopy(Matrix* dest, const Matrix* src) {
  if (!dest || !src) {
    return -1;
  }
  if ((dest->rows != src->rows) || (dest->cols != src->cols)) {
    fprintf(stderr, "Can't copy matrices of different sizes.\n");
    return -1;
  }
  memcpy(dest->data, src->data, slap_MatrixNumElements(dest) * sizeof(double));  // NOLINT
  return 0;
}

int slap_MatrixCopyFromArray(Matrix* mat, const double* data) {
  if (!mat) {
    return -1;
  }
  int len = slap_MatrixNumElements(mat);
  for (int i = 0; i < len; ++i) {
    mat->data[i] = data[i];
  }
  return 0;
}

int slap_MatrixCopyTranspose(Matrix* dest, Matrix* src) {
  if (!dest || !src) {
    return -1;
  }
  if ((dest->rows != src->cols) || (dest->cols != src->rows)) {
    fprintf(stderr,
            "Matrix sizes are not transposes of each other. Got (%d,%d) and (%d,%d).\n",
            dest->rows, dest->cols, src->rows, src->cols);
    return -1;
  }
  for (int i = 0; i < dest->rows; ++i) {
    for (int j = 0; j < dest->cols; ++j) {
      int dest_index = slap_MatrixGetLinearIndex(dest, i, j);
      int src_index = slap_MatrixGetLinearIndex(src, j, i);
      dest->data[dest_index] = src->data[src_index];
    }
  }
  return 0;
}

int slap_MatrixScaleByConst(Matrix* mat, double alpha) {
  if (!mat) {
    return -1;
  }
  for (int i = 0; i < slap_MatrixNumElements(mat); ++i) {
    mat->data[i] *= alpha;
  }
  return 0;
}

double slap_MatrixNormedDifference(const Matrix* A, const Matrix* B) {
  if (!A || !B) {
    return INFINITY;
  }
  if ((A->rows != B->rows) || (A->cols != B->cols)) {
    fprintf(stderr, "Can't compare matrices of different sizes. Got (%d,%d) and (%d,%d)\n",
            A->rows, A->cols, B->rows, B->cols);
    return INFINITY;
  }

  double diff = 0;
  for (int i = 0; i < slap_MatrixNumElements(A); ++i) {
    double d = A->data[i] - B->data[i];
    diff += d * d;
  }
  return sqrt(diff);
}

int slap_MatrixFlatten(Matrix* mat) {
  if (!mat) {
    return -1;
  }
  int size = slap_MatrixNumElements(mat);
  mat->rows = size;
  mat->cols = 1;
  return 0;
}

int slap_MatrixFlattenToRow(Matrix* mat) {
  if (!mat) {
    return -1;
  }
  int size = slap_MatrixNumElements(mat);
  mat->rows = 1;
  mat->cols = size;
  return 0;
}

int slap_PrintMatrix(const Matrix* mat) {
  if (!mat) {
    return -1;
  }
  for (int row = 0; row < mat->rows; ++row) {
    for (int col = 0; col < mat->cols; ++col) {
      printf("% 6.*g ", PRECISION, *slap_MatrixGetElementConst(mat, row, col));
    }
    printf("\n");
  }
  return 0;
}

int slap_PrintRowVector(const Matrix* mat) {
  if (!mat) {
    return -1;
  }
  printf("[ ");
  for (int i = 0; i < slap_MatrixNumElements(mat); ++i) {
    printf("% 6.*g ", PRECISION, mat->data[i]);
  }
  printf("]\n");
  return 0;
}

int slap_SetMatrixSize(Matrix* mat, int rows, int cols) {
  if (!mat) {
    return -1;
  }
  if (rows < 1 || cols < 1) {
    printf("ERROR: rows and columns must be positive integers.\n");
    return -1;
  }
  mat->rows = rows;
  mat->cols = cols;
  return 0;
}

int slap_MatrixSetIdentity(Matrix* mat, double val) {
  // TODO: verify square
  // TODO: check pointer
  slap_MatrixSetConst(mat, 0.0);
  for (int i = 0; i < mat->rows; ++i) {
    slap_MatrixSetElement(mat, i, i, val);
  }
  return 0;
}

int slap_MatrixSetDiagonal(Matrix* mat, const double* diag) {
  // TODO: check size (cols >= rows)
  for (int i = 0; i < mat->rows; ++i) {
    slap_MatrixSetElement(mat, i, i, diag[i]);
  }
  return 0;
}
