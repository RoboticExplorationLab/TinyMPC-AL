#pragma once
#include "matrix.h"

typedef struct {
  int rows;
  int cols;
  double* data;
  int stride_rows;
  int stride_cols;
} SubMatrix;

/**
 * @brief Create a submatrix of a matrix, or a "view" into a block of the 
 *        original matrix data
 * 
 * Note that the data is always assumed to be "owned" by a parent matrix, since the
 * data isn't assumed to be a contiguous block of memory. Hence, there are no "New" or 
 * "Free" methods for this type.
 * 
 * @param row  row of top-left corner of the submatrix 
 * @param col  column of the top-left corner of the submatrix 
 * @param rows number of rows in the view
 * @param cols number of columns in the view
 * @param A    original matrix
 * @return
 */
SubMatrix slap_SubMatrixFromMatrix(int row, int col, int rows, int cols, Matrix* A);

int slap_SubMatrixNumElements(const SubMatrix* mat);

int slap_SubMatrixCopyFromMatrix(SubMatrix* dest, const Matrix* src);

int slap_SubMatrixCopyToMatrix(Matrix* dest, const SubMatrix* src);

const double* slap_SubMatrixGetElementConst(const SubMatrix* mat, int row, int col);

double* slap_SubMatrixGetElement(SubMatrix* mat, int row, int col);

int slap_SubMatrixSetElement(SubMatrix* mat, int row, int col, double val);

int slap_SubMatrixSetConst(SubMatrix* mat, double val);

int slap_SubMatrixSetIdentity(SubMatrix* mat, double val);

int slap_SubMatrixCopyWithScaling(SubMatrix* dest, const Matrix* src, double val);
