#include "unconstrained_lqr.h"

void tiny_AddStageCost(double* cost, const tiny_ProblemData prob,
                       const Matrix x, const Matrix u, const int k) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[k], -1);
  *cost += 0.5 * slap_QuadraticForm(dx, prob.Q, dx);
  Matrix du = slap_MatrixFromArray(prob.ninputs, 1, dx_data);
  slap_MatrixAddition(du, u, prob.U_ref[k], -1);
  *cost += 0.5 * slap_QuadraticForm(du, prob.R, du);
}

void tiny_AddTerminalCost(double* cost, const tiny_ProblemData prob,
                          const Matrix x) {
  double dx_data[prob.nstates];
  Matrix dx = slap_MatrixFromArray(prob.nstates, 1, dx_data);
  slap_MatrixAddition(dx, x, prob.X_ref[prob.nhorizon - 1], -1);
  *cost += 0.5 * slap_QuadraticForm(dx, prob.Qf, dx);
}

void tiny_ExpandStageCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                          Matrix* hes_el_uu, Matrix* grad_el_u,
                          const tiny_ProblemData prob, const int k) {
  slap_MatrixCopy(*hes_el_xx, prob.Q);
  slap_MatMulAdd(*grad_el_x, prob.Q, prob.X_ref[k], -1, 0);
  slap_MatrixCopy(*hes_el_uu, prob.R);
  slap_MatMulAdd(*grad_el_u, prob.R, prob.U_ref[k], -1, 0);
}

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob) {
  slap_MatrixCopy(*hes_el_xx, prob.Qf);
  slap_MatMulAdd(*grad_el_x, prob.Qf, prob.X_ref[prob.nhorizon-1], -1, 0);
}

// Riccati recursion for LTI without constraints
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtiModel model,
                                         Matrix Q_temp) {
  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);

  Matrix Qxx = slap_CreateSubMatrix(Q_temp, 0, 0, n, n);
  Matrix Qxu = slap_CreateSubMatrix(Q_temp, 0, n, n, m);
  Matrix Qux = slap_CreateSubMatrix(Q_temp, n, 0, m, n);
  Matrix Quu = slap_CreateSubMatrix(Q_temp, n, n, m, m);
  Matrix Qx = slap_CreateSubMatrix(Q_temp, 0, n + m, n, 1);
  Matrix Qu = slap_CreateSubMatrix(Q_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);
  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, k);
    // State Gradient: Qx = q + A'(P*f + p)
    slap_MatrixCopy(prob->p[k], prob->p[k + 1]);
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
    slap_MatMulAdd(Quu, slap_Transpose(model.B), model.B, solver.regu, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Qux, slap_Transpose(model.B), prob->P[k], 1,
                   0);  // Qux = B'P*A

    // Calculate Gains
    slap_MatrixCopy(Quu_temp, Quu);
    slap_Cholesky(Quu_temp);
    slap_MatrixCopy(prob->K[k], Qux);
    slap_MatrixCopy(prob->d[k], Qu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Quu\Qu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Quu\Qux

    // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
    slap_MatrixCopy(prob->P[k], Qxx);                            // P = Qxx
    slap_MatMulAdd(Qxu, slap_Transpose(prob->K[k]), Quu, 1, 0);  // Qxu = K'Quu
    slap_MatMulAdd(prob->P[k], Qxu, prob->K[k], 1, 1);           // P += K'Quu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Qux, -2,
                   1);  // P -= K'Qux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Qux), prob->K[k], -1,
    //                1);  // P -= Qux'K

    // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
    slap_MatrixCopy(prob->p[k], Qx);                    // p = Qx
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

// Roll out the closed-loop dynamics with K and d from backward pass,
// calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtiModel model) {
  int N = prob.nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    // delta_x and delta_u over previous X, U
    // Control input: u = uf - d - K*(x - xf)
    slap_MatrixCopy(U[k], prob.d[k]);              // u[k] = d[k]
    slap_MatMulAdd(U[k], prob.K[k], X[k], -1, -1);   // u[k] += K[k] * x[k]
    // Next state: x = A*x + B*u + f
    tiny_LtiDynamics(&X[k + 1], X[k], U[k], model);
  }
  return SLAP_NO_ERROR;
}

// Roll out the closed-loop dynamics with K and d from backward pass,
// provided function to calculate Ak, Bk. Calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPassLtv(
    Matrix* X, Matrix* U, const tiny_ProblemData prob, 
    const tiny_LtvModel model) {
  int N = prob.nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    // Control input: u = - d - K*x
    slap_MatrixCopy(U[k], prob.d[k]);              // u[k] = -d[k]
    slap_MatMulAdd(U[k], prob.K[k], X[k], -1, -1);   // u[k] -= K[k] * x[k]
    // Next state: x = A*x + B*u + f
    tiny_LtvDynamics(&X[k + 1], X[k], U[k], model, k);
  }
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_BackwardPassLtv(
    tiny_ProblemData* prob, const tiny_Solver solver,
    const tiny_LtvModel model, Matrix Q_temp) {

  int N = prob->nhorizon;
  int n = prob->nstates;
  int m = prob->ninputs;
  tiny_ExpandTerminalCost(&(prob->P[N - 1]), &(prob->p[N - 1]), *prob);

  Matrix Qxx = slap_CreateSubMatrix(Q_temp, 0, 0, n, n);
  Matrix Qxu = slap_CreateSubMatrix(Q_temp, 0, n, n, m);
  Matrix Qux = slap_CreateSubMatrix(Q_temp, n, 0, m, n);
  Matrix Quu = slap_CreateSubMatrix(Q_temp, n, n, m, m);
  Matrix Qx = slap_CreateSubMatrix(Q_temp, 0, n + m, n, 1);
  Matrix Qu = slap_CreateSubMatrix(Q_temp, n, n + m, m, 1);

  // Hijack the first part of P[N] to use for Cholesky decomposition
  // NOTE: Assumes m <= n
  Matrix Quu_temp = slap_Reshape(prob->P[N - 1], m, m);
  for (int k = N - 2; k >= 0; --k) {
    // Stage cost expansion
    tiny_ExpandStageCost(&Qxx, &Qx, &Quu, &Qu, *prob, k);
    // State Gradient: Qx = q + A'(P*f + p)
    slap_MatrixCopy(prob->p[k], prob->p[k + 1]);
    slap_MatMulAdd(prob->p[k], prob->P[k + 1], model.f[k], 1,
                   1);  // p[k] = P[k+1]*f + p[k+1]
    slap_MatMulAdd(Qx, slap_Transpose(model.A[k]), prob->p[k], 1, 1);

    // Control Gradient: Qu = r + B'(P*f + p)
    slap_MatMulAdd(Qu, slap_Transpose(model.B[k]), prob->p[k], 1, 1);

    // State Hessian: Qxx = Q + A'P*A
    slap_MatMulAdd(prob->P[k], prob->P[k + 1], model.A[k], 1,
                   0);  // P[k] = P[k+1]*A
    slap_MatMulAdd(Qxx, slap_Transpose(model.A[k]), prob->P[k], 1,
                   1);  // Qxx = Q + A'P*A

    // Control Hessian Quu = R + B'P*B
    slap_MatMulAdd(Qxu, prob->P[k + 1], model.B[k], 1, 0);       // Qxu = P * B
    slap_MatMulAdd(Quu, slap_Transpose(model.B[k]), Qxu, 1, 1);  // Quu = R + B'P*B
    slap_MatMulAdd(Quu, slap_Transpose(model.B[k]), model.B[k], solver.regu, 1);
    // Hessian Cross-Term
    slap_MatMulAdd(Qux, slap_Transpose(model.B[k]), prob->P[k], 1,
                   0);  // Qux = B'P*A

    // Calculate Gains
    slap_MatrixCopy(Quu_temp, Quu);
    slap_Cholesky(Quu_temp);
    slap_MatrixCopy(prob->K[k], Qux);
    slap_MatrixCopy(prob->d[k], Qu);
    slap_CholeskySolve(Quu_temp, prob->d[k]);  // d = Quu\Qu
    slap_CholeskySolve(Quu_temp, prob->K[k]);  // K = Quu\Qux

    // Cost-to-Go Hessian: P = Qxx + K'Quu*K - K'Qux - Qux'K
    slap_MatrixCopy(prob->P[k], Qxx);                            // P = Qxx
    slap_MatMulAdd(Qxu, slap_Transpose(prob->K[k]), Quu, 1, 0);  // Qxu = K'Quu
    slap_MatMulAdd(prob->P[k], Qxu, prob->K[k], 1, 1);           // P += K'Quu*K
    slap_MatMulAdd(prob->P[k], slap_Transpose(prob->K[k]), Qux, -2,
                   1);  // P -= K'Qux
    // slap_MatMulAdd(prob->P[k], slap_Transpose(Qux), prob->K[k], -1,
    //                1);  // P -= Qux'K

    // Cost-to-Go Gradient: p = Qx + K'Quu*d - K'Qu - Qux'd
    slap_MatrixCopy(prob->p[k], Qx);                    // p = Qx
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

void tiny_LtiDynamics(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtiModel model) {
  slap_MatrixCopy(*xn, model.f);                          
  slap_MatMulAdd(*xn, model.A, x, 1, 1);      // x[k+1] += A * x[k]
  slap_MatMulAdd(*xn, model.B, u, 1, 1);      // x[k+1] += B * u[k]
}

void tiny_LtvDynamics(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtvModel model, const int k) {
  slap_MatrixCopy(*xn, model.f[k]);                          
  slap_MatMulAdd(*xn, model.A[k], x, 1, 1);      // x[k+1] += A * x[k]
  slap_MatMulAdd(*xn, model.B[k], u, 1, 1);      // x[k+1] += B * u[k]
}