//
// Created by Brian Jackson on 12/19/22.
// Copyright (c) 2022 Robotic Exploration Lab. All rights reserved.
//

#include "cholesky.h"

#include <math.h>

enum slap_ErrorCode slap_Cholesky(Matrix A) {
  SLAP_ASSERT_VALID(A, SLAP_INVALID_MATRIX, "Cholesky: matrix invalid");
  int n = slap_MinDim(A);
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < j; ++k) {
      for (int i = j; i < n; ++i) {
        double* Aij = slap_GetElement(A, i, j);
        double Aik = *slap_GetElement(A, i, k);
        double Ajk = *slap_GetElement(A, j, k);
        *Aij -= Aik * Ajk;
      }
    }
    double Ajj = *slap_GetElement(A, j, j);
    if (Ajj <= 0) {
      return SLAP_CHOLESKY_FAIL;
    }
    double ajj = sqrt(Ajj);

    for (int i = j; i < n; ++i) {
      double* Aij = slap_GetElement(A, i, j);
      *Aij /= ajj;
    }
  }
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode slap_TriSolve(Matrix L, Matrix b) {
  SLAP_ASSERT_VALID(L, SLAP_INVALID_MATRIX, "LowerTriBackSub: L matrix invalid");
  SLAP_ASSERT_VALID(b, SLAP_INVALID_MATRIX, "LowerTriBackSub: b matrix invalid");
  SLAP_ASSERT(slap_NumCols(L) == slap_NumRows(b), SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
              SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
              "LowerTriBackSub: L has %d columns but b has %d rows", slap_NumCols(L),
              slap_NumRows(b));
  int n = b.rows;
  int m = b.cols;
  bool tL = L.mattype == slap_TRIANGULAR_UPPER || slap_IsTransposed(L);

  for (int j_ = 0; j_ < n; ++j_) {
    int j = tL ? n - j_ - 1 : j_;
    for (int k = 0; k < m; ++k) {
      double* xjk = slap_GetElement(b, j, k);
      double Ljj = *slap_GetElement(L, j, j);
      *xjk /= Ljj;

      for (int i_ = j_ + 1; i_ < n; ++i_) {
        int i = tL ? i_ - (j_ + 1) : i_;
        double* xik = slap_GetElement(b, i, k);
        double Lij = *slap_GetElement(L, i, j);
        *xik -= Lij * (*xjk);
      }
    }
  }
  return SLAP_NO_ERROR;
}
enum slap_ErrorCode slap_CholeskySolve(const Matrix A, Matrix b) {
  // NOTE: Validity checks are done by the sub-methods
  enum slap_ErrorCode err;
  err = slap_TriSolve(A, b);
  if (err != SLAP_NO_ERROR) return err;
  err = slap_TriSolve(slap_Transpose(A), b);
  if (err != SLAP_NO_ERROR) return err;
  return SLAP_NO_ERROR;
}
