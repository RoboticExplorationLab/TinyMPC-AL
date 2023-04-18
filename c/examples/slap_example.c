#include "slap/slap.h"

int main() {

    sfloat A_data[4] = {1, 2, 3, 4};
    Matrix A;
    A = slap_MatrixFromArray(2, 2, A_data);
    // slap_PrintMatrix(A);

    sfloat B_data[4*2] = {2.0000000000000164,      0.4999999999999965,     2.000000000000016,
                          0.49999999999999645,     2.000000000000015,      0.4999999999999964,
                          2.0000000000000147,      0.49999999999999634};
    Matrix B[2];
    B[0] = slap_MatrixFromArray(2, 2, &B_data[0]);
    B[1] = slap_MatrixFromArray(2, 2, &B_data[4]);
    slap_PrintMatrix(B[0]);
    slap_PrintMatrix(B[1]);

    return 0;
}