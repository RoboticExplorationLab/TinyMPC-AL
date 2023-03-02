//
// Created by Brian Edward Jackson on 1/31/23.
// Check Lec 10 on DDP

#include "riccati.h"

enum slap_ErrorCode tiny_Riccati_LTI(
    int N, const Matrix A, const Matrix B,
    const Matrix Q, const Matrix R, const Matrix q,
    const Matrix r, Matrix* K, Matrix* d, Matrix* P,
    Matrix* p, Matrix S_temp) 
{
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
    // State Gradient: Sx = q + A' * p
    slap_MatrixCopy(p[k], p[k+1]);
    slap_MatrixCopy(Sx, q);                             // Sx = q
    slap_MatMulAdd(Sx, slap_Transpose(A), p[k], 1, 1);  // Sx += A' * p

    // Control Gradient: Su = r + B' * p
    slap_MatrixCopy(Su, r);                               // Su = r
    slap_MatMulAdd(Su, slap_Transpose(B), p[k], 1, 1);    // Su += B' * p

    // State Hessian: Sxx = Q + A'P*A
    slap_MatMulAdd(P[k], P[k + 1], A, 1, 0);  // P[k] = P * A (temp in P)
    slap_MatrixCopy(Sxx, Q);                  // Sxx = Q
    slap_MatMulAdd(Sxx, slap_Transpose(A), P[k], 1, 1);  // Sxx = Q + A'P*A

    // Control Hessian Suu = R + B'P*B
    slap_MatMulAdd(Sxu, P[k + 1], B, 1, 0);             // Sxu = P * B
    slap_MatrixCopy(Suu, R);                            // Suu = R
    slap_MatMulAdd(Suu, slap_Transpose(B), Sxu, 1, 1);  // Suu = R + B'P*B

    // Hessian Cross-Term
    slap_MatMulAdd(Sux, slap_Transpose(B), P[k], 1, 0);  // Sux = B'P*A

    // Calculate Gains
    slap_MatrixCopy(Quu_temp, Suu);
    slap_Cholesky(Quu_temp);
    slap_MatrixCopy(K[k], Sux);
    slap_MatrixCopy(d[k], Su);
    slap_CholeskySolve(Quu_temp, d[k]);  // d = Suu\Su
    slap_CholeskySolve(Quu_temp, K[k]);  // K = Suu\Sux

    // Cost-to-Go Hessian: P = Sxx + K'Suu*K - K'Sux - Sux'K
    slap_MatrixCopy(P[k], Sxx);                              // P = Sxx
    slap_MatMulAdd(Sxu, slap_Transpose(K[k]), Suu, 1, 0);    // Sxu = K'Suu
    slap_MatMulAdd(P[k], Sxu, K[k], 1, 1);                   // P += K'Suu*K
    slap_MatMulAdd(P[k], slap_Transpose(K[k]), Sux, -1, 1);   // P -= K'Sux
    slap_MatMulAdd(P[k], slap_Transpose(Sux), K[k], -1, 1);   // P -= Sux'K

    // Cost-to-Go Gradient: p = Sx + K'Suu*d - K'Su - Sux'd
    slap_MatrixCopy(p[k], Sx);                               // p = Sx
    slap_MatMulAdd(p[k], Sxu, d[k], 1, 1);                   // p += K'Suu*d
    slap_MatMulAdd(p[k], slap_Transpose(K[k]), Su, -1, 1);    // p -= K'Su
    slap_MatMulAdd(p[k], slap_Transpose(Sux), d[k], -1, 1);   // p -= Sux'd
  }
  slap_MatrixCopy(P[N], Q);  // Replace P[N] since we used it for Quu_temp
  // slap_FreeMatrix(temp);
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_Riccati_LTV(
    int N, const Matrix* A, const Matrix* B,
    const Matrix Q, const Matrix R, const Matrix q,
    const Matrix r, Matrix* K, Matrix* d, Matrix* P,
    Matrix* p, Matrix S_temp) 
{
  // Copy terminal cost-to-go
  int k = N;
  slap_MatrixCopy(P[k], Q);
  slap_MatrixCopy(p[k], q);

  int n = slap_NumCols(A[k]);
  int m = slap_NumCols(B[k]);

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
    // State Gradient: Sx = q + A' * p
    slap_MatrixCopy(p[k], p[k+1]);
    slap_MatrixCopy(Sx, q);                             // Sx = q
    slap_MatMulAdd(Sx, slap_Transpose(A[k]), p[k], 1, 1);  // Sx += A' * p

    // Control Gradient: Su = r + B' * p
    slap_MatrixCopy(Su, r);                               // Su = r
    slap_MatMulAdd(Su, slap_Transpose(B[k]), p[k], 1, 1);    // Su += B' * p

    // State Hessian: Sxx = Q + A'P*A
    slap_MatMulAdd(P[k], P[k+1], A[k], 1, 0);                   // P[k] = P * A (temp in P)
    slap_MatrixCopy(Sxx, Q);                                 // Sxx = Q
    slap_MatMulAdd(Sxx, slap_Transpose(A[k]), P[k], 1, 1);      // Sxx = Q + A'P*A

    // Control Hessian Suu = R + B'P*B
    slap_MatMulAdd(Sxu, P[k+1], B[k], 1, 0);                    // Sxu = P * B
    slap_MatrixCopy(Suu, R);                                 // Suu = R
    slap_MatMulAdd(Suu, slap_Transpose(B[k]), Sxu, 1, 1);       // Suu = R + B'P*B

    // Hessian Cross-Term
    slap_MatMulAdd(Sux, slap_Transpose(B[k]), P[k], 1, 0);       // Sux = B'P*A

    // Calculate Gains
    slap_MatrixCopy(Quu_temp, Suu);
    slap_Cholesky(Quu_temp);
    slap_MatrixCopy(K[k], Sux);
    slap_MatrixCopy(d[k], Su);
    slap_CholeskySolve(Quu_temp, d[k]);  // d = Suu\Su
    slap_CholeskySolve(Quu_temp, K[k]);  // K = Suu\Sux

    // Cost-to-Go Hessian: P = Sxx + K'Suu*K - K'Sux - Sux'K
    slap_MatrixCopy(P[k], Sxx);                              // P = Sxx
    slap_MatMulAdd(Sxu, slap_Transpose(K[k]), Suu, 1, 0);    // Sxu = K'Suu
    slap_MatMulAdd(P[k], Sxu, K[k], 1, 1);                   // P += K'Suu*K
    slap_MatMulAdd(P[k], slap_Transpose(K[k]), Sux, -1, 1);   // P -= K'Sux
    slap_MatMulAdd(P[k], slap_Transpose(Sux), K[k], -1, 1);   // P -= Sux'K

    // Cost-to-Go Gradient: p = Sx + K'Suu*d - K'Su - Sux'd
    slap_MatrixCopy(p[k], Sx);                               // p = Sx
    slap_MatMulAdd(p[k], Sxu, d[k], 1, 1);                   // p += K'Suu*d
    slap_MatMulAdd(p[k], slap_Transpose(K[k]), Su, -1, 1);    // p -= K'Su
    slap_MatMulAdd(p[k], slap_Transpose(Sux), d[k], -1, 1);   // p -= Sux'd

  }
  slap_MatrixCopy(P[N], Q);  // Replace P[N] since we used it for Quu_temp
  // slap_FreeMatrix(temp);
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_Riccati_LTVf(
    int N, Matrix A, Matrix B,
    void (*get_jacobians)(Matrix, Matrix, const Matrix, const Matrix), 
    const Matrix Q, const Matrix R, const Matrix q, const Matrix r, 
    Matrix* K, Matrix* d, Matrix* P, Matrix* p, 
    const Matrix* x_ref, const Matrix* u_ref, Matrix S_temp)
{
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
    // Get Jacobians of dynamics
    (*get_jacobians)(A, B, x_ref[k], u_ref[k]);
    // State Gradient: Sx = q + A'(P*f + p)
    slap_MatrixCopy(p[k], p[k+1]);
    // slap_PrintMatrix(p[k]);
    slap_MatrixCopy(Sx, q);                             // Sx = q
    
    slap_MatMulAdd(Sx, slap_Transpose(A), p[k], 1, 1);  // Sx += A' * p
    
    // Control Gradient: Su = r + B' * p
    slap_MatrixCopy(Su, r);                               // Su = r
    slap_MatMulAdd(Su, slap_Transpose(B), p[k], 1, 1);    // Su += B' * p

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

    // Cost-to-Go Hessian: P = Sxx + K'Suu*K - K'Sux - Sux'K
    slap_MatrixCopy(P[k], Sxx);                              // P = Sxx
    slap_MatMulAdd(Sxu, slap_Transpose(K[k]), Suu, 1, 0);    // Sxu = K'Suu
    slap_MatMulAdd(P[k], Sxu, K[k], 1, 1);                   // P += K'Suu*K
    slap_MatMulAdd(P[k], slap_Transpose(K[k]), Sux, -1, 1);   // P -= K'Sux
    slap_MatMulAdd(P[k], slap_Transpose(Sux), K[k], -1, 1);   // P -= Sux'K

    // Cost-to-Go Gradient: p = Sx + K'Suu*d - K'Su - Sux'd
    slap_MatrixCopy(p[k], Sx);                               // p = Sx
    slap_MatMulAdd(p[k], Sxu, d[k], 1, 1);                   // p += K'Suu*d
    slap_MatMulAdd(p[k], slap_Transpose(K[k]), Su, -1, 1);    // p -= K'Su
    slap_MatMulAdd(p[k], slap_Transpose(Sux), d[k], -1, 1);   // p -= Sux'd
  }
  slap_MatrixCopy(P[N], Q);  // Replace P[N] since we used it for Quu_temp
  // slap_FreeMatrix(temp);
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_RiccatiForwardPass_LTI(int N, const Matrix A, const Matrix B,
                                            const Matrix x0, const Matrix xf, const Matrix uf, 
                                            const Matrix* K, const Matrix* d,
                                            const Matrix* P, const Matrix* p, Matrix* x,
                                            Matrix* u, Matrix* y) 
{
  // Initial condition
  slap_MatrixCopy(x[0], x0);

  for (int k = 0; k <= N; ++k) {
    // Calculate primal variables (via closed-loop forward simulation)
    if (k < N) {
      // Control input: u = uf - d - K*(x - xf) 
      slap_MatrixAddition(u[k], uf, d[k], -1);    // u[k] = uf + d[k]
      slap_MatMulAdd(u[k], K[k], x[k], -1, 1);   // u[k] -= K[k] * x[k]
      slap_MatMulAdd(u[k], K[k], xf, 1, 1);   // u[k] += K[k] * xf
      // Next state: x = A*x + B*u
      slap_MatMulAdd(x[k + 1], A, x[k], 1, 0);  // x[k+1] = A * x[k]
      slap_MatMulAdd(x[k + 1], B, u[k], 1, 1);  // x[k+1] += B * u[k]
    }

    // Calculate dual variables
    if (y) {
      slap_MatrixCopy(y[k], p[k]);             // y[k] = p[k];
      slap_MatMulAdd(y[k], P[k], x[k], 1, 1);  // y[k] += P[k] x[k]
    }
  }
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_RiccatiForwardPass_LTV(int N, const Matrix* A, const Matrix* B,
                                            const Matrix x0, const Matrix xf, const Matrix uf, 
                                            const Matrix* K, const Matrix* d,
                                            const Matrix* P, const Matrix* p, Matrix* x,
                                            Matrix* u, Matrix* y) 
{
  // Initial condition
  slap_MatrixCopy(x[0], x0);

  for (int k = 0; k <= N; ++k) {
    // Calculate primal variables (via closed-loop forward simulation)
    if (k < N) {
      // Control input: u = uf - d - K*(x - xf) 
      slap_MatrixAddition(u[k], uf, d[k], -1);    // u[k] = uf + d[k]
      slap_MatMulAdd(u[k], K[k], x[k], -1, 1);   // u[k] -= K[k] * x[k]
      slap_MatMulAdd(u[k], K[k], xf, 1, 1);   // u[k] += K[k] * xf
      // Next state: x = A*x + B*u
      slap_MatMulAdd(x[k + 1], A[k], x[k], 1, 0);  // x[k+1] = A * x[k]
      slap_MatMulAdd(x[k + 1], B[k], u[k], 1, 1);  // x[k+1] += B * u[k]
    }

    // Calculate dual variables
    if (y) {
      slap_MatrixCopy(y[k], p[k]);              // y[k] = p[k];
      slap_MatMulAdd(y[k], P[k], x[k], 1, 1);   // y[k] += P[k] x[k]
    }
  }
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_RiccatiForwardPass_LTVf(
        int N, Matrix A, Matrix B,
        void (*get_jacobians)(Matrix, Matrix, const Matrix, const Matrix), 
        const Matrix x0, const Matrix* xref, const Matrix* uref, 
        const Matrix* K, const Matrix* d, const Matrix* P, const Matrix* p, 
        Matrix* x, Matrix* u, Matrix* y)
{
  // Initial condition
  slap_MatrixCopy(x[0], x0);

  for (int k = 0; k <= N; ++k) {
    // Calculate primal variables (via closed-loop forward simulation)
    if (k < N) {
      // Get Jacobians of dynamics
      (*get_jacobians)(A, B, xref[k], uref[k]);
      // Control input: u = uf - d - K*(x - xf) 
      slap_MatrixAddition(u[k], uref[k], d[k], -1);    // u[k] = uf - d[k]
      slap_MatMulAdd(u[k], K[k], x[k], -1, 1);   // u[k] -= K[k] * x[k]
      slap_MatMulAdd(u[k], K[k], xref[k], 1, 1);   // u[k] += K[k] * xf
      // Next state: x = A*x + B*u + f
      tiny_Dynamics_RK4_Raw(x[k+1].data, x[k].data, u[k].data);
      // slap_MatMulAdd(x[k + 1], A, x[k], 1, 0);  // x[k+1] = A * x[k]
      // slap_MatMulAdd(x[k + 1], B, u[k], 1, 1);  // x[k+1] += B * u[k]
    }

    // Calculate dual variables
    if (y) {
      slap_MatrixCopy(y[k], p[k]);              // y[k] = p[k];
      slap_MatMulAdd(y[k], P[k], x[k], 1, 1);   // y[k] += P[k] x[k]
    }
  }
  return SLAP_NO_ERROR;  
}

// enum slap_ErrorCode tiny_LQR_LTI(int N, const Matrix A, const Matrix B, 
//                                  const Matrix Q, const Matrix R, const Matrix q,
//                                  const Matrix r, Matrix* K, Matrix* d, Matrix* P,
//                                  Matrix* p, Matrix S_temp)
// {
//   // Copy terminal cost-to-go
//   int k = N;
//   slap_MatrixCopy(P[N], Q);
//   slap_MatrixCopy(p[N], q);

//   int n = slap_NumCols(A);
//   int m = slap_NumCols(B);

//   // Create a new local matrix variable
//   Matrix temp_nn = slap_CreateSubMatrix(S_temp, 0, 0, n, n);
//   Matrix temp_nm = slap_CreateSubMatrix(S_temp, 0, n, n, m);
//   Matrix temp_mn = slap_CreateSubMatrix(S_temp, n, 0, m, n);
//   Matrix temp_mm = slap_CreateSubMatrix(S_temp, n, n, m, m);
//   Matrix temp_nn2 = slap_Reshape(P[N], n, n);  //hijack P[N]
//   // Matrix temp_nn2 = slap_NewMatrix(n, n);
//   for (--k; k >= 0; --k) 
//   {
//     // K = (R + B' * P * B) \ (B' * P * A)
//     slap_MatMulAdd(temp_nm, P[k+1], B, 1, 0);  //temp_nm = P * B
//     slap_MatrixCopy(temp_mm, R);
//     slap_MatMulAdd(temp_mm, slap_Transpose(B), temp_nm, 1, 1); //R + B' * P * B
//     slap_MatMulAdd(temp_nn, P[k+1], A, 1, 0);  //temp_nn = P * A
//     slap_MatMulAdd(K[k], slap_Transpose(B), temp_nn, 1, 0); //B' * P * A
//     slap_Cholesky(temp_mm);
//     slap_CholeskySolve(temp_mm, K[k]);
    
//     // P = Q + A' * P * (A - B * K)
//     slap_MatrixCopy(temp_nn, A);
//     slap_MatMulAdd(temp_nn, B, K[k], -1, 1);  //A - B * K
//     if (k == N-1) //not affect using P[N] data
//       slap_MatMulAdd(temp_nn2, Q, temp_nn, 1, 0);  //P * (A - B * K)
//     else
//       slap_MatMulAdd(temp_nn2, P[k+1], temp_nn, 1, 0);  //P * (A - B * K)
//     slap_MatrixCopy(P[k], Q);
//     slap_MatMulAdd(P[k], slap_Transpose(A), temp_nn2, 1, 1);
//   }
//   // slap_FreeMatrix(temp_nn2);
//   return SLAP_NO_ERROR;
// }                      

// enum slap_ErrorCode tiny_LQR_LTV(int N, const Matrix* A, const Matrix* B, 
//                                  const Matrix Q, const Matrix R, const Matrix q,
//                                  const Matrix r, Matrix* K, Matrix* d, Matrix* P,
//                                  Matrix* p, Matrix S_temp)
// {
//   // Copy terminal cost-to-go
//   int k = N;
//   slap_MatrixCopy(P[N], Q);
//   slap_MatrixCopy(p[N], q);

//   int n = slap_NumCols(A[0]);
//   int m = slap_NumCols(B[0]);

//   // Create a new local matrix variable
//   Matrix temp_nn = slap_CreateSubMatrix(S_temp, 0, 0, n, n);
//   Matrix temp_nm = slap_CreateSubMatrix(S_temp, 0, n, n, m);
//   Matrix temp_mn = slap_CreateSubMatrix(S_temp, n, 0, m, n);
//   Matrix temp_mm = slap_CreateSubMatrix(S_temp, n, n, m, m);
//   Matrix temp_nn2 = slap_Reshape(P[N], n, n);  //hijack P[N]
//   // Matrix temp_nn2 = slap_NewMatrix(n, n);
//   for (--k; k >= 0; --k) 
//   {
//     // K = (R + B' * P * B) \ (B' * P * A)
//     slap_MatMulAdd(temp_nm, P[k+1], B[k], 1, 0);  //temp_nm = P * B
//     slap_MatrixCopy(temp_mm, R);
//     slap_MatMulAdd(temp_mm, slap_Transpose(B[k]), temp_nm, 1, 1); //R + B' * P * B
//     slap_MatMulAdd(temp_nn, P[k+1], A[k], 1, 0);  //temp_nn = P * A
//     slap_MatMulAdd(K[k], slap_Transpose(B[k]), temp_nn, 1, 0); //B' * P * A
//     slap_Cholesky(temp_mm);
//     slap_CholeskySolve(temp_mm, K[k]);
    
//     // P = Q + A' * P * (A - B * K)
//     slap_MatrixCopy(temp_nn, A[k]);
//     slap_MatMulAdd(temp_nn, B[k], K[k], -1, 1);  //A - B * K
//     if (k == N-1) //not affect using P[N] data
//       slap_MatMulAdd(temp_nn2, Q, temp_nn, 1, 0);  //P * (A - B * K)
//     else
//       slap_MatMulAdd(temp_nn2, P[k+1], temp_nn, 1, 0);  //P * (A - B * K)
//     slap_MatrixCopy(P[k], Q);
//     slap_MatMulAdd(P[k], slap_Transpose(A[k]), temp_nn2, 1, 1);
//   }
//   // slap_FreeMatrix(temp_nn2);
//   return SLAP_NO_ERROR;
// }

// enum slap_ErrorCode tiny_LQR_LTVf(
//         int N, Matrix A, Matrix B,
//         void (*get_jacobians)(Matrix, Matrix, const Matrix, const Matrix), 
//         const Matrix Q, const Matrix R, const Matrix q, const Matrix r, 
//         Matrix* K, Matrix* d, Matrix* P, Matrix* p, 
//         const Matrix* Un, const Matrix* Xn, Matrix S_temp)
// {
//   // Copy terminal cost-to-go
//   int k = N;
//   slap_MatrixCopy(P[N], Q);
//   slap_MatrixCopy(p[N], q);

//   int n = slap_NumCols(A);
//   int m = slap_NumCols(B);

//   // Create a new local matrix variable
//   Matrix temp_nn = slap_CreateSubMatrix(S_temp, 0, 0, n, n);
//   Matrix temp_nm = slap_CreateSubMatrix(S_temp, 0, n, n, m);
//   Matrix temp_mn = slap_CreateSubMatrix(S_temp, n, 0, m, n);
//   Matrix temp_mm = slap_CreateSubMatrix(S_temp, n, n, m, m);
//   Matrix temp_nn2 = slap_Reshape(P[N], n, n);  //hijack P[N]
//   // Matrix temp_nn2 = slap_NewMatrix(n, n);
//   for (--k; k >= 0; --k) 
//   {
//     (*get_jacobians)(A, B, Xn[k], Un[k]);
//     // K = (R + B' * P * B) \ (B' * P * A)
//     slap_MatMulAdd(temp_nm, P[k+1], B, 1, 0);  //temp_nm = P * B
//     slap_MatrixCopy(temp_mm, R);
//     slap_MatMulAdd(temp_mm, slap_Transpose(B), temp_nm, 1, 1); //R + B' * P * B
//     slap_MatMulAdd(temp_nn, P[k+1], A, 1, 0);  //temp_nn = P * A
//     slap_MatMulAdd(K[k], slap_Transpose(B), temp_nn, 1, 0); //B' * P * A
//     slap_Cholesky(temp_mm);
//     slap_CholeskySolve(temp_mm, K[k]);
    
//     // P = Q + A' * P * (A - B * K)
//     slap_MatrixCopy(temp_nn, A);
//     slap_MatMulAdd(temp_nn, B, K[k], -1, 1);  //A - B * K
//     if (k == N-1) //not affect using P[N] data
//       slap_MatMulAdd(temp_nn2, Q, temp_nn, 1, 0);  //P * (A - B * K)
//     else
//       slap_MatMulAdd(temp_nn2, P[k+1], temp_nn, 1, 0);  //P * (A - B * K)
//     slap_MatrixCopy(P[k], Q);
//     slap_MatMulAdd(P[k], slap_Transpose(A), temp_nn2, 1, 1);
//   }
//   // slap_FreeMatrix(temp_nn2);
//   return SLAP_NO_ERROR;
// }