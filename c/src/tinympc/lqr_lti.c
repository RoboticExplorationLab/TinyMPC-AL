#include "lqr_lti.h"

// Riccati recursion for LTI without constraints
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtiModel model,
                                         Matrix* Q_temp) {
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);

  Matrix Qxx = slap_CreateSubMatrix(*Q_temp, 0, 0, n, n);
  Matrix Qxu = slap_CreateSubMatrix(*Q_temp, 0, n, n, m);
  Matrix Qux = slap_CreateSubMatrix(*Q_temp, n, 0, m, n);
  Matrix Quu = slap_CreateSubMatrix(*Q_temp, n, n, m, m);
  Matrix Qx = slap_CreateSubMatrix(*Q_temp, 0, n + m, n, 1);
  Matrix Qu = slap_CreateSubMatrix(*Q_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);
  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, k);
    // State Gradient: Qx = q + A'(P*f + p)
    slap_Copy(prob->p[k], prob->p[k + 1]);
    slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f, 1,
                   1);  // p[k] = P[k+1]*f + p[k+1]
    slap_MatMulAdd(Qx, slap_Transpose(model.A), prob->p[k], 1, 1);

    // Control Gradient: Qu = r + B'(P*f + p)
    slap_MatMulAdd(Qu, slap_Transpose(model.B), prob->p[k], 1, 1);

    // State Hessian: Qxx = Q + A'P*A
    slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A, 1,
                   0);  // P[k] = P[k+1]*A
    slap_MatMulAdd(Qxx, slap_Transpose(model.A), prob->P[k], 1,
                   1);  // Qxx = Q + A'P*A

    // Control Hessian Quu = R + B'P*B
    slap_MatMulAdd(Qxu, prob->P[k + 1], model.B, 1, 0);       // Qxu = P * B
    slap_MatMulAdd(Quu, slap_Transpose(model.B), Qxu, 1, 1);  // Quu = R + B'P*B
    slap_MatMulAdd(Quu, slap_Transpose(model.B), model.B, solver.reg, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Qux, slap_Transpose(model.B), prob->P[k], 1,
                   0);  // Qux = B'P*A

    // Calculate Gains
    slap_Copy(Quu_temp, Quu);
    slap_Cholesky(Quu_temp);
    slap_Copy(prob->K[k], Qux);
    slap_Copy(prob->d[k], Qu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Quu\Qu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Quu\Qux

    // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
    slap_Copy(prob->P[k], Qxx);                                  // P = Qxx
    slap_MatMulAdd(Qxu, slap_Transpose(prob->K[k]), Quu, 1, 0);  // Qxu = K'Quu
    slap_MatMulAdd(prob->P[k], Qxu, prob->K[k], 1, 1);           // P += K'Quu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Qux, -2,
                   1);  // P -= K'Qux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Qux), prob->K[k], -1,
    //                1);  // P -= Qux'K

    // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
    slap_Copy(prob->p[k], Qx);                          // p = Qx
    slap_MatMulAdd(prob->p[k], Qxu, prob->d[k], 1, 1);  // p += K'Quu*d
    slap_MatMulAdd(prob->p[k], slap_Transpose(prob->K[k]), Qu, -1,
                   1);  // p -= K'Qu
    slap_MatMulAdd(prob->p[k], slap_Transpose(Qux), prob->d[k], -1,
                   1);  // p -= Qux'd
  }
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);
  // Replace P[N] since we used it for Quu_temp (need improving later)
  return SLAP_NO_ERROR;
}