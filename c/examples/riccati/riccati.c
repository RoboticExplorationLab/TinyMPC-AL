//
// Created by Brian Edward Jackson on 1/31/23.
//

#include "riccati.h"
enum slap_ErrorCode slap_Riccati_LTI(int N, const Matrix A, const Matrix B, const Matrix f,
                                     const Matrix Q, const Matrix R, const Matrix q,
                                     const Matrix r, Matrix* K, Matrix* d, Matrix* P,
                                     Matrix* p, Matrix S_temp) {
  // Copy terminal cost-to-go
  int k = N;
  slap_MatrixCopy(P[k], Q);
  slap_MatrixCopy(p[k], q);

  int n = slap_NumCols(A);
  int m = slap_NumCols(B);

  Matrix Sxx = slap_CreateSubMatrix(S_temp, 0, 0, n, n);
  Matrix Sxu = slap_CreateSubMatrix(S_temp, 0, n, n, m);
  Matrix Sux = slap_CreateSubMatrix(S_temp, n, 0, m, n);
  Matrix Suu = slap_CreateSubMatrix(S_temp, n, n, m, m);
  Matrix Sx = slap_CreateSubMatrix(S_temp, 0, n + m, n, 1);
  Matrix Su = slap_CreateSubMatrix(S_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(P[N], m, m);

  for (--k; k >= 0; --k) {
    // State Gradient: Sx = q + A' * (P * f + p)
    slap_MatrixCopy(p[k], p[k+1]);
    slap_MatMulAdd(p[k], P[k+1], f, 1, 1);              // p = P * f + p
    slap_MatrixCopy(Sx, q);                             // Sx = q
    slap_MatMulAdd(Sx, slap_Transpose(A), p[k], 1, 1);  // Sx = q + A' * (P*f + p)

    // Control Gradient: Su = r + B' * (P * f + p)
    slap_MatrixCopy(Su, r);                               // Su = r
    slap_MatMulAdd(Su, slap_Transpose(B), p[k], 1, 1);    // Su = r + B' * (P*f + p)

    // State Hessian: Sxx = Q + A'P*A
    slap_MatMulAdd(P[k], P[k+1], A, 1, 0);                   // P[k] = P * A (temp in P)
    slap_MatrixCopy(Sxx, Q);                                 // Sxx = Q
    slap_MatMulAdd(Sxx, slap_Transpose(A), P[k], 1, 1);      // Sxx = Q + A'P*A

    // Control Hessian Suu = R + B'P*B
    slap_MatMulAdd(Sxu, P[k+1], B, 1, 0);                    // Sxu = P * B
    slap_MatrixCopy(Suu, R);                                 // Suu = R
    slap_MatMulAdd(Suu, slap_Transpose(B), Sxu, 1, 1);       // Suu = R + B'P*B

    // Hessian Cross-Term
    slap_MatMulAdd(Sux, slap_Transpose(B), P[k], 1, 0);       // Sux = B'P*A

    // Calculate Gains
    slap_MatrixCopy(Quu_temp, Suu);
    slap_Cholesky(Quu_temp);
    slap_MatrixCopy(K[k], Sux);
    slap_MatrixCopy(d[k], Su);
    slap_CholeskySolve(Quu_temp, d[k]);  // d = Suu\Su
    slap_CholeskySolve(Quu_temp, K[k]);  // K = Suu\Sux
    slap_ScaleByConst(K[k], -1);
    slap_ScaleByConst(d[k], -1);

    // Cost-to-Go Hessian: P = Sxx + K'Suu*K + K'Sux + Sux'K
    slap_MatrixCopy(P[k], Sxx);                              // P = Sxx
    slap_MatMulAdd(Sxu, slap_Transpose(K[k]), Suu, 1, 0);    // Sxu = K'Suu
    slap_MatMulAdd(P[k], Sxu, K[k], 1, 1);                   // P += K'Suu*K
    slap_MatMulAdd(P[k], slap_Transpose(K[k]), Sux, 1, 1);   // P += K'Sux
    slap_MatMulAdd(P[k], slap_Transpose(Sux), K[k], 1, 1);   // P += Sux'K

    // Cost-to-Go Gradient: p = Sx + K'Suu*d + K'Su + Sux'd
    slap_MatrixCopy(p[k], Sx);                               // p = Sx
    slap_MatMulAdd(p[k], Sxu, d[k], 1, 1);                   // p += K'Suu*d
    slap_MatMulAdd(p[k], slap_Transpose(K[k]), Su, 1, 1);    // p += K'Su
    slap_MatMulAdd(p[k], slap_Transpose(Sux), d[k], 1, 1);   // p += Sux'd
  }
  slap_MatrixCopy(P[N], Q);  // Replace P[N] since we used it for Quu_temp

//  slap_FreeMatrix(temp);
  return SLAP_NO_ERROR;
}
enum slap_ErrorCode slap_RiccatiForwardPass_LTI(int N, const Matrix A, const Matrix B,
                                            const Matrix f, const Matrix x0,
                                            const Matrix* K, const Matrix* d,
                                            const Matrix* P, const Matrix* p, Matrix* x,
                                            Matrix* u, Matrix* y) {
  // Initial condition
  slap_MatrixCopy(x[0], x0);

  for (int k = 0; k <= N; ++k) {
    // Calculate primal variables (via closed-loop forward simulation)
    if (k < N) {
      slap_MatrixCopy(u[k], d[k]);              // u[k] = d[k]
      slap_MatMulAdd(u[k], K[k], x[k], 1, 1);   // u[k] += K[k] * x[k]
      slap_MatrixCopy(x[k + 1], f);             // x[k+1] = f
      slap_MatMulAdd(x[k + 1], A, x[k], 1, 1);  // x[k+1] += A * x[k]
      slap_MatMulAdd(x[k + 1], B, u[k], 1, 1);  // x[k+1] += B * u[k]
    }

    // Calculate dual variables
    if (y) {
      slap_MatrixCopy(y[k], p[k]);              // y[k] = p[k];
      slap_MatMulAdd(y[k], P[k], x[k], 1, 1);   // y[k] += P[k] x[k]
    }
  }
  return SLAP_NO_ERROR;
}
