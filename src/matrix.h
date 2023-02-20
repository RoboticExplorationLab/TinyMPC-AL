/**
 * @file matrix.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Matrix type and basic operations
 * @version 0.1
 * @date 2022-01-30
 *
 * @copyright Copyright (c) 2022
 *
 * @addtogroup LinearAlgebra
 * @{
 */
#pragma once

#include <stdbool.h>

/**
 * @brief Represents a matrix of double-precision data
 *
 * Simple wrapper around an arbitrary pointer to the underlying data.
 * The data is assumed to be stored in a contiguous block of memory.
 * The data is interpreted column-wise, such that `data[1]` is element `[1,0]` of the
 * matrix.
 *
 * ## Initialization
 * A Matrix can be initialized a few ways. The easiest is via `NewMatrix`:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Matrix mat = slap_NewMatrix(rows, cols);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * which allocates a new block of memory on the heap. It must be followed by a call to
 * FreeMatrix().
 *
 * If the data for the matrix is already stored in an array, the default brace initializer
 * can be used:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * double data[6] = {1,2,3,4,5,6};
 * Matrix mat = {2, 3, data};
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * ## Methods
 * The following methods are defined for the Matrix type:
 *
 * ### Initialization and deconstruction
 * - slap_NewMatrix()
 * - FreeMatrix()
 * - MatrixSetConst()
 * - MatrixScaleByConst()
 *
 * ### Indexing operations
 * - MatrixNumElements()
 * - MatrixGetLinearIndex()
 * - MatrixGetElement()
 * - MatrixGetElementTranspose()
 * - MatrixSetElement()
 *
 * ### Copying
 * - MatrixCopy()
 * - MatrixCopyTranspose()
 *
 * ### Reshaping
 * - MatrixFlatten()
 * - MatrixFlattenToRow()
 *
 * ### Utilities
 * - PrintMatrix()
 * - PrintRowVector()
 * - MatrixNormedDifference()
 */
typedef struct {
  int rows;
  int cols;
  double* data;
} Matrix;

/**
 * @brief Allocate a new matrix on the heap
 *
 * Data will not be initialized. Wrapper around a call to `malloc`.
 * Must be followed by a call to `FreeMatrix`.
 *
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 * @return A new matrix
 */
Matrix slap_NewMatrix(int rows, int cols);

/**
 * @brief Wraps existing data in a Matrix class
 * 
 * @param rows Number of rows in the matrix
 * @param cols Number of columns in the matrix
 * @param data Data for the matrix. Must not be NULL, and should have at least 
 *             rows * cols elements.
 * @return A new matrix
 */
Matrix slap_MatrixFromArray(int rows, int cols, double* data);

/**
 * @brief Allocate a new matrix on the heap, initialized with zeros
 *
 * Data will not be initialized. Wrapper around a call to `malloc`.
 * Must be followed by a call to `FreeMatrix`.
 *
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 * @return A new matrix
 */
Matrix slap_NewMatrixZeros(int rows, int cols);

/**
 * @brief Sets all of the elements in a matrix to a single value
 *
 * @param mat Matrix to be modified
 * @param val Value to which each element will be set
 * @return 0 if successful
 */
int slap_MatrixSetConst(Matrix* mat, double val);

/**
 * @brief Free the data for a matrix
 *
 * Note this does NOT attempt to free the matrix object itself, only the data
 * it wraps.
 *
 * @param mat
 * @post [mat.data](Matrix.data) will be `NULL`.
 * @return 0 if successful
 */
int slap_FreeMatrix(Matrix* mat);

/**
 * @brief Get the number of elements in a matrix, i.e. `m * n`.
 *
 * @param mat Any matrix
 * @return Number of elements in the matrix
 */
int slap_MatrixNumElements(const Matrix* mat);

/**
 * @brief Get the linear index for a given row and column in the matrix
 *
 * Converts a cartesian index of row and column into a linear index for accessing
 * an element of the underlying data.
 *
 * @param mat Matrix with nonzero size and initialized data
 * @param row Row index
 * @param col Column index
 * @return Linear index corresponding to `row` and `col`.
           Returns -1 for a bad input.
 */
int slap_MatrixGetLinearIndex(const Matrix* mat, int row, int col);

/**
 * @brief Get the element of a matrix or its transpose
 *
 * If @p istransposed is false, then this method acts just like MatrixGetElement().
 * Otherwise, it is equalivalent to flipping the @p row and @p col arguments to
 * MatrixGetElement().
 *
 * @param mat         Matrix with nonzero size and initialized data
 * @param row         Row index
 * @param col         Column index
 * @param istranposed Are the indicies for the transpose of A?
 * @return            Pointer to the data at the given element.
 */
double* slap_MatrixGetElementTranspose(Matrix* mat, int row, int col,
                                       bool istranposed);

const double* slap_MatrixGetElementTransposeConst(const Matrix* mat, int row, int col,
                                                  bool istranposed);
/**
 * @brief The a matrix element to a given value
 *
 * @param mat Matrix with nonzero size and initialized data
 * @param row Row index
 * @param col Column index
 * @param val Value to which the element should be set
 * @return    0 if successful
 */
int slap_MatrixSetElement(Matrix* mat, int row, int col, double val);

/**
 * @brief Get the element of a matrix given row, column indices
 *
 * @param mat Matrix of nonzero size
 * @param row Row index
 * @param col Column index
 * @return A pointer to the element of the matrix. NULL for invalid input.
 */
double* slap_MatrixGetElement(Matrix* mat, int row, int col);

const double* slap_MatrixGetElementConst(const Matrix* mat, int row, int col);

/**
 * @brief Copy a matrix to another matrix, transposed
 *
 * @param dest a matrix of size (m,n)
 * @param src a matrix of size (n,m)
 * @return 0 if successful
 */
int slap_MatrixCopyTranspose(Matrix* dest, Matrix* src);

/**
 * @brief Copy a matrix to another matrix
 *
 * @param dest a matrix of size (m,n)
 * @param src a matrix of size (n,m)
 * @return 0 if successful
 */
int slap_MatrixCopy(Matrix* dest, const Matrix* src);

/**
 * @brief Copy the data from an array into the matrix, columnwise.
 *
 * @param mat  Matrix with nonzero size
 * @param data Data to be copied into the array. Must have length of at least mat.rows *
 * mat.cols.
 * @return 0 if successful
 */
int slap_MatrixCopyFromArray(Matrix* mat, const double* data);

/**
 * @brief Scale a matrix by a constant factor
 *
 * @param mat Fully initialized matrix of non-zero size. Values will be modified.
 * @param alpha scalar by which to multiply the matrix
 * @return 0 if successsful
 */
int slap_MatrixScaleByConst(Matrix* mat, double alpha);

/**
 * @brief Return the normed difference between 2 matrices of the same size
 *
 * Returns \f$ \sqrt{\sum_{i=0}^{m-1} \sum_{j=0}^{n-1} (A_{ij} - B_{ij})^2 } \f$
 *
 * @param A A matrix of dimension (m,n)
 * @param B A matrix of dimension (m,n)
 * @return
 */
double slap_MatrixNormedDifference(const Matrix* A, const Matrix* B);

/**
 * @brief Flatten a 2D matrix to a column vector
 *
 * Changes the row and column data so that the matrix is now a column vector. The
 * underlying data is unchanged.
 *
 * @param mat Matrix to be flattened.
 * @return 0 if successful
 */
int slap_MatrixFlatten(Matrix* mat);

/**
 * @brief Flatten a 2D matrix to a row vector
 *
 * Changes the row and column data so that the matrix is now a row vector. The
 * underlying data is unchanged.
 *
 * @param mat Matrix to be flattened
 * @return 0 if successful
 */
int slap_MatrixFlattenToRow(Matrix* mat);

/**
 * @brief Print the elements of a matrix to stdout
 *
 * Precision of the printing can be controlled by the global variable PRECISION.
 *
 * @param mat Matrix to be printed
 * @return 0 if successful
 */
int slap_PrintMatrix(const Matrix* mat);

/**
 * @brief Print the entire matrix as a row vector
 *
 * Same result as calling PrintMatrix() after a call to MatrixFlattenToRow().
 *
 * @param mat Matrix to be printed
 * @return 0 if successful
 */
int slap_PrintRowVector(const Matrix* mat);

/**
 * @brief Set the dimensions of the matrix
 *
 * Note that this does not change the underlying data, only it's interpretation.
 *
 * @param mat  Matrix
 * @param rows New number of rows
 * @param cols New number of columns
 * @return 0 if successful
 */
int slap_SetMatrixSize(Matrix* mat, int rows, int cols);

/**
 * @brief Set the diagonal elements of the matrix to val, and the rest to zeros.
 *
 * @param mat Square matrix
 * @param val Value for the diagonal elements
 * @return
 */
int slap_MatrixSetIdentity(Matrix* mat, double val);

/**
 * @brief Set the matrix diagonal from an array
 *
 * Doesn't touch any of the off-diagonal elements.
 *
 * @param mat Matrix (nrows >= ncols)
 * @param diag Array of length `nrows`.
 * @return
 */
int slap_MatrixSetDiagonal(Matrix* mat, const double* diag);


/**@}*/
