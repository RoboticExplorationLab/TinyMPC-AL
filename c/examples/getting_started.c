#include <math.h>
#include <stdlib.h>

#include "slap/slap.h"

double myfun(double x) { return 2 * x * sin(x); }

int main(void) {
  puts("Welcome to Getting Started with the slap Library!");

  /////////////////////////////////////////////
  // Basic Operations
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~~~ BASIC OPS ~~~~~~~~~~~~~~\n");

  // Create a matrix with stack-allocated memory
  double data_A[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix A = slap_MatrixFromArray(3, 4, data_A);

  // Get the sizes of the matrix
  int n_rows = slap_NumRows(A);
  int n_cols = slap_NumCols(A);
  int n_el = slap_NumElements(A);
  printf("Size of A = (%d,%d) with %d elements\n", n_rows, n_cols, n_el);

  // Get a pointer to an element
  double *pa = slap_GetElement(A, 0, 0);
  double a = *slap_GetElement(A, 0, 1);
  printf("A[0,0] = %0.3g\n", *pa);
  printf("A[0,1] = %0.3g\n", a);

  // Set an element via pointer
  *pa = -20;
  printf("A[0,0] = %0.3g (after set)\n", *slap_GetElement(A, 0, 0));

  // Set an element directly
  slap_SetElement(A, 1, 2, 25.5);
  printf("A[1,2] = %0.3g (after set)\n", *slap_GetElement(A, 1, 2));

  /////////////////////////////////////////////
  // Unary Operations
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~~~ UNARY OPS ~~~~~~~~~~~~~~\n");

  // Printing
  printf("\nA:\n");
  slap_PrintMatrix(A);

  // Set the matrix to a constant
  const double val = 1.2;
  slap_SetConst(A, val);
  printf("\nA: (set const)\n");
  slap_PrintMatrix(A);

  // Scale by constant
  slap_ScaleByConst(A, -1);
  printf("\nA: (Scale by const)\n");
  slap_PrintMatrix(A);

  // Set Identity
  slap_SetIdentity(A, 2.1);
  printf("\nA: (Set to identity)\n");
  slap_PrintMatrix(A);

  // Set Diagonal
  double diag[3] = {1, 2, 3};
  slap_SetDiagonal(A, diag, 3);
  printf("\nA: (Set Diagonal)\n");
  slap_PrintMatrix(A);

  // Set first n elements of diagonal
  diag[0] = 10;
  diag[1] = 11;
  slap_SetDiagonal(A, diag, 2);
  printf("\nA: (Set Partial Diagonal)\n");
  slap_PrintMatrix(A);

  // Set Range
  slap_SetRange(A, 1, slap_NumElements(A));
  printf("\nA: (Set Range)\n");
  slap_PrintMatrix(A);

  slap_SetRange(A, 0, 1);
  printf("\nA: (Set Range 0 to 1)\n");
  slap_PrintMatrix(A);

  /////////////////////////////////////////////
  // Transformations
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~ TRANSFORMATIONS ~~~~~~~~~~~\n");

  // Transpose
  Matrix At = slap_Transpose(A);
  printf("\nA transpose:\n");
  slap_PrintMatrix(At);
  bool tA = slap_IsTransposed(At);
  printf("At is transposed? %d\n", tA);

  // Flatten
  Matrix vec_a = slap_Flatten(A);
  printf("\nA Flatten\n");
  slap_PrintMatrix(vec_a);

  // Reshape
  Matrix A2 = slap_Reshape(A, 2, 6);
  printf("\nA2\n");
  slap_PrintMatrix(A2);

  // Modify reshape (note that it changes the original)
  slap_SetElement(A2, 0, 1, 100);
  printf("\nA (after edit via reshape)\n");
  slap_PrintMatrix(A);

  // Reshape to smaller
  // note this is 1st 4 elements, not top left corner
  Matrix A_resize = slap_Reshape(A, 2, 2);
  printf("\nA resize\n");
  slap_PrintMatrix(A_resize);

  /////////////////////////////////////////////
  // Copying
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~~~ COPYING ~~~~~~~~~~~~~~~\n");

  // Matrix with heap-allocated memory
  double *data_B = (double *)malloc(n_el * sizeof(double));
  Matrix B = slap_MatrixFromArray(4, 3, data_B);

  // Copy from transposed array
  slap_MatrixCopy(B, At);

  printf("\nCopy from A to B:\n");
  printf("A:\n");
  slap_PrintMatrix(A);
  printf("B:\n");
  slap_PrintMatrix(B);

  // Copy from array
  double data_C[4] = {-1, 2, -3, 4};
  slap_MatrixCopyFromArray(
      A_resize, data_C);  // note we're copying to a reshaped version of A

  printf("\nA (after array copy):\n");
  slap_PrintMatrix(A);

  /////////////////////////////////////////////
  // Sub-Arrays
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~~ SUB-ARRAYS ~~~~~~~~~~~~~~\n");

  // Get view of top-left 2x2 corner
  Matrix A_sub = slap_CreateSubMatrix(A, 0, 0, 2, 2);
  printf("\nTop-left corner of A\n");
  slap_PrintMatrix(A_sub);
  printf("A is Dense?     %d\n", slap_IsDense(A));
  printf("A sub is Dense? %d\n", slap_IsDense(A_sub));

  // Copy to Sub-matrix
  data_C[3] = -50;
  slap_MatrixCopyFromArray(A_sub, data_C);
  printf("\nA (after copying to sub-matrix)\n");
  slap_PrintMatrix(A);

  // Get middle elements
  A_sub = slap_CreateSubMatrix(A, 1, 1, 1, 2);
  printf("\nMiddle of A\n");
  slap_PrintMatrix(A_sub);

  /////////////////////////////////////////////
  // Vector Operations
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~~ VECTOR OPS ~~~~~~~~~~~~~~\n");

  // Create Matrices (vectors) on the heap
  // Note the calls to slap_FreeMatrix at the bottom of this function
  Matrix x = slap_NewMatrix(3, 1);
  Matrix y = slap_NewMatrix(3, 1);
  Matrix z = slap_NewMatrixZeros(3, 1);

  // Set some values
  slap_SetRange(x, 1, 3);
  slap_SetConst(y, 1.5);

  // Inner product
  double dot_xy = slap_InnerProduct(x, y);
  printf("Inner product = %0.3g\n", dot_xy);

  // Cross product
  slap_CrossProduct(z, x, y);
  printf("Cross product: ");
  slap_PrintMatrix(slap_Transpose(z));

  // Norms (only for dense)
  printf("Two Norm Squared: %0.2g\n", slap_NormTwoSquared(z));
  printf("Two Norm:         %0.2g\n", slap_NormTwo(z));
  printf("One Norm:         %0.2g\n", slap_NormOne(z));
  printf("Inf Norm:         %0.2g\n", slap_NormInf(z));

  // Max/Min
  printf("Max: %0.2g\n", slap_Max(z));
  printf("Min: %0.2g\n", slap_Min(z));

  // Argmax/Argmin
  double z_max, z_min;
  MatrixIterator it_max = slap_ArgMax(z, &z_max);
  MatrixIterator it_min = slap_ArgMin(z, &z_min);
  printf("ArgMax: %0.2g at linear index, %d, Cartesian index (%d,%d)\n", z_max,
         it_max.k, it_max.i, it_max.j);
  printf("ArgMin: %0.2g at linear index, %d, Cartesian index (%d,%d)\n", z_max,
         it_min.k, it_min.i, it_min.j);

  // Outer Product
  Matrix Z = slap_NewMatrix(3, 3);
  slap_OuterProduct(Z, x, y);
  printf("\nOuter Product:\n");
  slap_PrintMatrix(Z);

  /////////////////////////////////////////////
  // Linear Algebra
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~~ LINEAR ALG ~~~~~~~~~~~~~~\n");

  // Reshape A to be square
  A = slap_Reshape(A, 3, 3);

  // Create output C matrix
  Matrix C = slap_NewMatrix(3, 4);
  slap_SetConst(C, 1);

  // Matrix Multiplication (C = beta * C + alpha * A * B)
  double alpha = 1.0;
  double beta = 0.5;
  slap_MatMulAdd(C, A, slap_Transpose(B), alpha, beta);
  printf("\nC (MatMulAdd)\n");
  slap_PrintMatrix(C);

  // Linear Solve w/ Cholesky
  double data_A2[9];                       // data for PD matrix
  double data_b[3] = {10.1, -11.2, 12.3};  // rhs vector
  A2 = slap_MatrixFromArray(3, 3, data_A2);
  Matrix b = slap_MatrixFromArray(3, 1, data_b);

  slap_MatMulAtB(A2, A,
                 A);  // A2 = A'A. NOTE: Use with caution. Please see docs.
  slap_AddIdentity(A2, 0.01);  // Ensure A2 > 0
  slap_MatrixCopy(A, A2);      // Save A2 back to A for later

  enum slap_ErrorCode err;
  err = slap_Cholesky(A2);  // Perform Cholesky decomposition
  if (err == SLAP_CHOLESKY_FAIL) {
    printf("Matrix is not Positive-Definite!\n");
  }

  slap_MatrixCopy(x, b);      // copy rhs to x
  slap_CholeskySolve(A2, x);  // solve for x

  slap_MatMulAB(y, A, x);  // Calculate y = A * x (y should equal b)
  printf("\nCholesky Solve (these should be equal):\n");
  printf("b = ");
  slap_PrintMatrix(slap_Transpose(b));
  printf("y = ");
  slap_PrintMatrix(slap_Transpose(y));

  /////////////////////////////////////////////
  // Advanced Ops
  /////////////////////////////////////////////

  printf("\n~~~~~~~~~~~~~ ADVANCED OPS ~~~~~~~~~~~~~\n");

  // Efficient Iteration (including Sub-Arrays)
  Matrix B_sub = slap_CreateSubMatrix(B, 1, 1, 3, 2);
  printf("\nB\n");
  slap_PrintMatrix(B);

  printf("\nB_sub\n");
  slap_PrintMatrix(B_sub);

  printf("\nIteration over Sub-Array:\n");
  for (MatrixIterator it = slap_Iterator(B_sub); !slap_IsFinished(&it);
       slap_Step(&it)) {
    int mem_index = it.index;
    int lin_index = it.k;
    int row_index = it.i;
    int col_index = it.j;
    double B_val = B_sub.data[mem_index];
    printf("B_sub[%d,%d] = % 4.2f at linear index %d and memory index %d\n",
           row_index, col_index, B_val, lin_index, mem_index);
  }

  // Function Mapping
  printf("\nB (before map)\n");
  slap_PrintMatrix(B);

  slap_Map(B_sub, myfun);

  printf("\nB (after map)\n");
  slap_PrintMatrix(B);

  /////////////////////////////////////////////
  // Memory Cleanup (DON'T FORGET THIS)
  /////////////////////////////////////////////

  free(data_B);
  slap_FreeMatrix(x);
  slap_FreeMatrix(y);
  slap_FreeMatrix(z);
  slap_FreeMatrix(Z);
  slap_FreeMatrix(C);
  return EXIT_SUCCESS;
}
