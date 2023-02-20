#include "submatrix.h"

SubMatrix slap_SubMatrixFromMatrix(int row, int col, int rows, int cols, Matrix* A) {
  // TODO: do size check 
  // TODO: do pointer check
  double* data_new = slap_MatrixGetElement(A, row, col);
  SubMatrix mat = {rows, cols, data_new, 1, A->rows};
  return mat;
}

int slap_SubMatrixNumElements(const SubMatrix* mat) {
  // TODO: do pointer check
  return mat->rows * mat->cols;
}

int slap_SubMatrixCopyFromMatrix(SubMatrix* dest, const Matrix* src) {
  // TODO: do size check
  // TODO: do pointer check
  for (int j = 0; j < dest->cols; ++j) {
    for (int i = 0; i < dest->rows; ++i) {
      double* el_dest = slap_SubMatrixGetElement(dest, i, j);
      const double* el_src = slap_MatrixGetElementConst(src, i, j);
      *el_dest = *el_src;
    }
  }
  return 0;
}

int slap_SubMatrixCopyToMatrix(Matrix* dest, const SubMatrix* src) {
  // TODO: do size check
  // TODO: do pointer check
  for (int j = 0; j < dest->cols; ++j) {
    for (int i = 0; i < dest->rows; ++i) {
      double* el_dest = slap_MatrixGetElement(dest, i, j);
      const double* el_src = slap_SubMatrixGetElementConst(src, i, j);
      *el_dest = *el_src;
    }
  }
  return 0;
}

const double* slap_SubMatrixGetElementConst(const SubMatrix* mat, int row, int col) {
  int index = row * mat->stride_rows + col * mat->stride_cols;
  return mat->data + index;
}

double* slap_SubMatrixGetElement(SubMatrix* mat, int row, int col) {
  int index = row * mat->stride_rows + col * mat->stride_cols;
  return mat->data + index;
}

int slap_SubMatrixSetElement(SubMatrix* mat, int row, int col, double val) {
  // TODO: check size and pointer
  *slap_SubMatrixGetElement(mat, row, col) = val;
  return 0;
}

int slap_SubMatrixSetConst(SubMatrix* mat, double val) {
  for (int j = 0; j < mat->cols; ++j) {
    for (int i = 0; i < mat->rows; ++i) {
      slap_SubMatrixSetElement(mat, i, j, val);
    }
  }
  return 0;
}

int slap_SubMatrixSetIdentity(SubMatrix* mat, double val) {
  // TODO: verify square
  // TODO: check pointer
  slap_SubMatrixSetConst(mat, 0.0);
  for (int i = 0; i < mat->rows; ++i) {
    slap_SubMatrixSetElement(mat, i, i, val);
  }
  return 0;
}
int slap_SubMatrixCopyWithScaling(SubMatrix* dest, const Matrix* src, double val) {
  // TODO: check sizes
  for (int j = 0; j < dest->cols; ++j) {
    for (int i = 0; i < dest->rows; ++i) {
      double* el_dest = slap_SubMatrixGetElement(dest, i, j);
      const double* el_src = slap_MatrixGetElementConst(src, i, j);
      *el_dest = *el_src * val;
    }
  }
  return 0;
}
