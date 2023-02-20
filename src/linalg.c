#include "linalg.h"

#include "math.h"
#include "matrix.h"
#include "stdio.h"

int slap_MatrixAddition(Matrix* C, const Matrix* A, const Matrix* B, double alpha) {
  for (int i = 0; i < slap_MatrixNumElements(C); ++i) {
    C->data[i] = A->data[i] + alpha * B->data[i];
  }
  return 0;
}

int slap_MatrixScale(Matrix* A, double alpha) {
  for (int i = 0; i < slap_MatrixNumElements(A); ++i) {
    A->data[i] *= alpha;
  }
  return 0;
}

int slap_MatrixMultiply(Matrix* C, const Matrix* A, const Matrix* B, bool tA, bool tB,
                        double alpha, double beta) {
  int n;
  int m;
  if (tA) {
    n = A->cols;
    m = A->rows;
  } else {
    n = A->rows;
    m = A->cols;
  }
  int p = tB ? B->rows : B->cols;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      double* Cij = slap_MatrixGetElement(C, i, j);
      *Cij *= beta;
      for (int k = 0; k < m; ++k) {
        double Aik = *slap_MatrixGetElementTransposeConst(A, i, k, tA);
        double Bkj = *slap_MatrixGetElementTransposeConst(B, k, j, tB);
        *Cij += alpha * Aik * Bkj;
      }
    }
  }
  return 0;
}

int slap_SymmetricMatrixMultiply(Matrix* Asym, Matrix* B, Matrix* C, double alpha,
                                 double beta) {
  int n;
  int m;
  bool tA = false;
  bool tB = false;
  if (tA) {
    n = Asym->cols;
    m = Asym->rows;
  } else {
    n = Asym->rows;
    m = Asym->cols;
  }
  int p = tB ? B->rows : B->cols;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      double* Cij = slap_MatrixGetElement(C, i, j);
      *Cij *= beta;
      for (int k = 0; k < m; ++k) {
        int row = i;
        int col = k;
        if (i < k) {
          row = k;
          col = i;
        }
        double Aik = *slap_MatrixGetElement(Asym, row, col);
        double Bkj = *slap_MatrixGetElement(B, k, j);
        *Cij += alpha * Aik * Bkj;
      }
    }
  }
  return 0;
  return 0;
}

int slap_AddIdentity(Matrix* A, double alpha) {
  int n = A->rows;
  for (int i = 0; i < n; ++i) {
    double* Aii = slap_MatrixGetElement(A, i, i);
    *Aii += alpha;
  }
  return 0;
}

int slap_CholeskyFactorize(Matrix* A) {
  int n = A->rows;
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < j; ++k) {
      for (int i = j; i < n; ++i) {
        double* Aij = slap_MatrixGetElement(A, i, j);
        double Aik = *slap_MatrixGetElement(A, i, k);
        double Ajk = *slap_MatrixGetElement(A, j, k);
        *Aij -= Aik * Ajk;
      }
    }
    double Ajj = *slap_MatrixGetElement(A, j, j);
    if (Ajj <= 0) {
      return slap_kCholeskyFail;
    }
    double ajj = sqrt(Ajj);

    for (int i = j; i < n; ++i) {
      double* Aij = slap_MatrixGetElement(A, i, j);
      *Aij /= ajj;
    }
  }
  return slap_kCholeskySuccess;
}

int slap_LowerTriBackSub(Matrix* L, Matrix* b, bool istransposed) {
  int n = b->rows;
  int m = b->cols;
  for (int j_ = 0; j_ < n; ++j_) {
    int j = istransposed ? n - j_ - 1 : j_;
    for (int k = 0; k < m; ++k) {
      double* xjk = slap_MatrixGetElement(b, j, k);
      double Ljj = *slap_MatrixGetElement(L, j, j);
      *xjk /= Ljj;

      for (int i_ = j_ + 1; i_ < n; ++i_) {
        int i = istransposed ? i_ - (j_ + 1) : i_;
        double* xik = slap_MatrixGetElement(b, i, k);
        double Lij = *slap_MatrixGetElementTranspose(L, i, j, istransposed);
        *xik -= Lij * (*xjk);
      }
    }
  }
  return 0;
}

int slap_CholeskySolve(Matrix* L, Matrix* b) {
  slap_LowerTriBackSub(L, b, 0);
  slap_LowerTriBackSub(L, b, 1);
  return 0;
}

double slap_TwoNorm(const Matrix* M) {
  if (!M) {
    return -1;
  }
  double norm = 0.0;
  for (int i = 0; i < slap_MatrixNumElements(M); ++i) {
    double x = M->data[i];
    norm += x * x;
  }
  return sqrt(norm);
}

double slap_OneNorm(const Matrix* M) {
  if (!M) {
    return -1;
  }
  double norm = 0.0;
  for (int i = 0; i < slap_MatrixNumElements(M); ++i) {
    double x = M->data[i];
    norm += fabs(x);
  }
  return sqrt(norm);
}

double slap_DotProduct(const Matrix* x, const Matrix* y) {
  if ((x->rows != y->rows) || (x->cols != 1) || (y->cols != 1)) {
    return NAN;
  }
  double out = 0.0;
  for (int i = 0; i < x->rows; ++i) {
    double xi = x->data[i];
    double yi = y->data[i];
    out += xi * yi;
  }
  return out;
}

double slap_QuadraticForm(const Matrix* x, const Matrix* A, const Matrix* y) {
  if ((x->rows != A->rows) || (y->rows != A->cols) || (x->cols != 1) || (y->cols != 1)) {
    return NAN;
  }
  double out = 0.0;
  for (int i = 0; i < x->rows; ++i) {
    for (int j = 0; j < y->rows; ++j) {
      double xi = x->data[i];
      double yj = y->data[j];
      double Aij = *slap_MatrixGetElementConst(A, i, j);
      out += xi * Aij * yj;
    }
  }
  return out;
}
double slap_DiagonalMultiplyLeft(Matrix* B, const Matrix* D, const Matrix* A) {
  // TODO: check sizes
  for (int j = 0; j < B->cols; ++j) {
    for (int i = 0; i < B->rows; ++i) {
      // TODO: do a single loop and use i = k % rows?
      double aij = *slap_MatrixGetElementConst(A, i, j);
      double di = D->data[i];
      slap_MatrixSetElement(B, i, j, di * aij);
    }
  }
  return 0;
}
double slap_DiagonalMultiplyRight(Matrix* B, const Matrix* A, const Matrix* D) {
  // TODO: check sizes
  for (int j = 0; j < B->cols; ++j) {
    double dj = D->data[j];
    for (int i = 0; i < B->rows; ++i) {
      // TODO: do a single loop and use i = k / cols?
      double aij = *slap_MatrixGetElementConst(A, i, j);
      slap_MatrixSetElement(B, i, j, dj * aij);
    }
  }
  return 0;
}
