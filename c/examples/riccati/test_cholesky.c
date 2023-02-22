#include <stdlib.h>
#include <math.h>
#include "slap/slap.h"

int main()
{
    // Linear Solve w/ Cholesky
    double data_A2[2] = {1., 0., 0., 1.};                       // data for PD matrix
    double data_b[3] = {10.1, -11.2, 12.3};  // rhs vector
    double data_A1[1] = {0.1};
    Matrix A2 = slap_MatrixFromArray(2, 2, data_A2);
    Matrix A1 = slap_MatrixFromArray(1, 1, data_A2);
    Matrix b = slap_MatrixFromArray(3, 1, data_b);
    slap_PrintMatrix(A1);
    slap_Cholesky(A1);
    slap_PrintMatrix(A1);
    slap_PrintMatrix(A2);
    slap_Cholesky(A2);
    slap_PrintMatrix(A2);
}
