#include "augmented_lagrangian_lqr.h"

const tiny_LinearDiscreteModel kDefaultLinearDiscreteModel = {
  .nstates = 0,
  .ninputs = 0,
  .dt = 0.0,
  .A = kNullMat,
  .B = kNullMat,
  .f = kNullMat,
  .x0 = kNullMat,
};

const tiny_KnotPoint kDefaultKnotPoint = {
  .x = kNullMat,
  .u = kNullMat,
  .t = 0.0,
  .dt = 0.0,
};

const tiny_Solver kDefaultSolver = {
  .regu = 0.0,
  .input_duals = NULL,
  .state_duals = NULL,
  .goal_dual = kNullMat,
  .penalty_min = 0.0,
  .penalty_max = 0.0,
  .penalty_mul = 0.0,
};

const tiny_ProblemData kDefaultProblemData = {
  .nstates = 0,
  .ninputs = 0,
  .nhorizon = 0,
  .ncstr_states = 0,
  .ncstr_inputs = 0,
  .ncstr_goal = 0,
  .Q = kNullMat,
  .R = kNullMat,
  .q = kNullMat,
  .r = kNullMat,
  .Qf = kNullMat,
  .u_max = kNullMat,
  .u_min = kNullMat,
  .x_max = kNullMat,
  .x_min = kNullMat,
  .X_ref = NULL,
  .U_ref = NULL,
  .dt = 0.0,
  .x0 = kNullMat,
  .K = NULL,
  .d = NULL,
  .P = NULL,
  .p = NULL,  
};

enum slap_ErrorCode tiny_BackwardPass(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    const Matrix* X, const Matrix* U, const tiny_Solver solver, Matrix S_temp) {
  // Copy terminal cost-to-go
  int k = prob.nhorizon;
  slap_MatrixCopy(P[k], prob.Qf);
  slap_MatrixCopy(p[k], prob.qf);
  int n = slap_NumCols(A);
  int m = slap_NumCols(B);

  return SLAP_NO_ERROR;
}

// Roll out the closed-loop dynamics with K and d from backward pass, 
// calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPass(
    const tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    Matrix* X, Matrix* U) {
  int N = prob.nhorizon;
  double alpha = 1.0;
  double u_data[prob.ninputs];
  double x_data[prob.nstates];
  Matrix Un = slap_MatrixFromArray(prob.ninputs, 1, u_data);
  Matrix Xn = slap_MatrixFromArray(prob.nstates, 1, x_data);
  slap_MatrixCopy(Xn, prob.x0);
  for (int k = 0; k < N - 1; ++k) {
    // delta_x and delta_u over previous X, U
    // Control input: u = uf - d - K*(x - xf) 
    slap_MatrixAddition(Un, U[k], prob.d[k], -1);    // u[k] = uf - d[k]
    slap_MatMulAdd(Un, prob.K[k], X[k], -1, 1);   // u[k] -= K[k] * x[k]
    slap_MatMulAdd(Un, prob.K[k], Xn, 1, 1);   // u[k] += K[k] * xf
    slap_MatrixCopy(U[k], Un);
    // Next state: x = A*x + B*u + f
    slap_MatrixCopy(Xn, X[k + 1]);
    slap_MatMulAdd(X[k + 1], model.A, X[k], 1, 0);  // x[k+1] = A * x[k]
    slap_MatMulAdd(X[k + 1], model.B, U[k], 1, 1);  // x[k+1] += B * u[k]
    slap_MatrixAddition(X[k + 1], X[k + 1], model.f, 1);  // x[k+1] += f
  }
  return SLAP_NO_ERROR; 
}

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model,
    Matrix* X, Matrix* U, const int verbose) {
  return SLAP_NO_ERROR;
}

// [u-p.u_max;-u + p.u_min]
void tiny_IneqInputs(
    Matrix ineq, const tiny_ProblemData prob, const Matrix u) {
  slap_SetConst(ineq, 0);  //clear before processing
  Matrix upper_half = slap_CreateSubMatrix(ineq, 0, 0, prob.ninputs, 1);
  Matrix lower_half = slap_CreateSubMatrix(ineq, prob.ninputs, 0, 
                                          prob.ninputs, 1);
  slap_MatrixAddition(upper_half, u, prob.u_max, -1);
  slap_MatrixAddition(lower_half, prob.u_min, u, -1);
}

void tiny_IneqInputsJacobian(
    Matrix ineq_jac, const tiny_ProblemData prob, const Matrix u) {
  slap_SetConst(ineq_jac, 0);  //clear before processing
  Matrix upper_half = slap_CreateSubMatrix(ineq_jac, 0, 0, 
                                          prob.ninputs, prob.ninputs);
  Matrix lower_half = slap_CreateSubMatrix(ineq_jac, prob.ninputs, 0, 
                                          prob.ninputs, prob.ninputs);  
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
}
  
void tiny_ActiveIneqMask(
    Matrix mask, const Matrix input_dual, const Matrix ineq) {    
  slap_SetConst(mask, 0);  //clear before processing
  for (int i = 0; i < ineq.rows; ++i) {
    // When variables are on the boundary or violating constraints
    bool active = input_dual.data[i] > 0 || ineq.data[i] > 0;
    slap_SetElement(mask, i, i, active); 
  }
}

