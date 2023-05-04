#include "constraint_linear.h"

enum tiny_ErrorCode tiny_SetInputBound(tiny_Workspace* work, sfloat* Ac_data, sfloat* bc_data) {
  int n = work->data->model->ninputs;
  work->stgs->en_cstr_inputs = EN_CSTR_INPUTS;
  work->data->Acu = slap_MatrixFromArray(2 * n, n, Ac_data);
  Matrix upper_half = slap_CreateSubMatrix(work->data->Acu, 0, 0, n, n);
  Matrix lower_half = slap_CreateSubMatrix(work->data->Acu, n, 0, n, n);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  work->data->bcu = slap_MatrixFromArray(2 * n, 1, bc_data);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetStateBound(tiny_Workspace* work, sfloat* Ac_data, sfloat* bc_data) {
  int n = work->data->model->nstates;
  work->stgs->en_cstr_states = EN_CSTR_STATES;
  work->data->Acx = slap_MatrixFromArray(2 * n, n, Ac_data);
  Matrix upper_half = slap_CreateSubMatrix(work->data->Acx, 0, 0, n, n);
  Matrix lower_half = slap_CreateSubMatrix(work->data->Acx, n, 0, n, n);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  work->data->bcx = slap_MatrixFromArray(2 * n, 1, bc_data);
  return TINY_NO_ERROR;
}

// // [u-p.u_max;-u + p.u_min]
// void tiny_EvalInputConstraint(Matrix* cu, const tiny_ProblemData prob,
//                      const Matrix u) {
//   slap_SetConst(*cu, 0);  // clear before processing
//   Matrix upper_half = slap_CreateSubMatrix(*cu, 0, 0, prob.ninputs,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cu, prob.ninputs, 0, prob.ninputs, 1);
//   slap_MatrixAddition(upper_half, u, prob.u_max, -1);
//   slap_MatrixAddition(lower_half, prob.u_min, u, -1);
// }

// Input constraints: A*u - b <= 0
enum tiny_ErrorCode tiny_EvalInputConstraint(tiny_Workspace* work, const int k) {
  // slap_SetConst(*cu, 0);  // clear before processing
  slap_MatMulAdd(work->cu, work->data->Acu, work->soln->U[k], 1.0, 0.0);
  slap_MatrixAddition(work->cu, work->cu, work->data->bcu, -1.0);
  return TINY_NO_ERROR;
}

// // [u_max, -u_min]
// void tiny_EvalInputConstraintOffset(Matrix* cu, const tiny_ProblemData prob) {
//   Matrix upper_half = slap_CreateSubMatrix(*cu, 0, 0, prob.ninputs,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cu, prob.ninputs, 0, prob.ninputs, 1);
//   slap_Copy(upper_half, prob.u_max);
//   slap_Copy(lower_half, prob.u_min);
//   slap_ScaleByConst(lower_half, -1);
// }

// void tiny_EvalInputConstraintJacobian(Matrix* ineq_jac, const tiny_ProblemData prob) {
//   slap_SetConst(*ineq_jac, 0);  // clear before processing
//   Matrix upper_half =
//       slap_CreateSubMatrix(*ineq_jac, 0, 0, prob.ninputs, prob.ninputs);
//   Matrix lower_half = slap_CreateSubMatrix(*ineq_jac, prob.ninputs, 0,
//                                            prob.ninputs, prob.ninputs);
//   slap_SetIdentity(upper_half, 1);
//   slap_SetIdentity(lower_half, -1);
// }

// // [x-p.x_max;-x + p.x_min]
// void tiny_EvalStateConstraint(Matrix* cx, const tiny_ProblemData prob,
//                      const Matrix x) {
//   slap_SetConst(*cx, 0);  // clear before processing
//   Matrix upper_half = slap_CreateSubMatrix(*cx, 0, 0, prob.nstates,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cx, prob.nstates, 0, prob.nstates, 1);
//   slap_MatrixAddition(upper_half, x, prob.x_max, -1);
//   slap_MatrixAddition(lower_half, prob.x_min, x, -1);
// }

// State constraints: A*x - b <= 0
enum tiny_ErrorCode tiny_EvalStateConstraint(tiny_Workspace* work, const int k) {
  // slap_SetConst(*cx, 0);  // clear before processing
  slap_MatMulAdd(work->cx, work->data->Acx, work->soln->X[k], 1.0, 0.0);
  slap_MatrixAddition(work->cx, work->cx, work->data->bcx, -1.0);
  return TINY_NO_ERROR;
}

// // [x_max, -x_min]
// void tiny_EvalStateConstraintOffset(Matrix* cx, const tiny_ProblemData prob) {
//   Matrix upper_half = slap_CreateSubMatrix(*cx, 0, 0, prob.nstates,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cx, prob.nstates, 0, prob.nstates, 1);
//   slap_Copy(upper_half, prob.x_max);
//   slap_Copy(lower_half, prob.x_min);
//   slap_ScaleByConst(lower_half, -1);
// }

// void tiny_EvalStateConstraintJacobian(Matrix* ineq_jac, const tiny_ProblemData prob) {
//   slap_SetConst(*ineq_jac, 0);  // clear before processing
//   Matrix upper_half =
//       slap_CreateSubMatrix(*ineq_jac, 0, 0, prob.nstates, prob.nstates);
//   Matrix lower_half = slap_CreateSubMatrix(*ineq_jac, prob.nstates, 0,
//                                            prob.nstates, prob.nstates);
//   slap_SetIdentity(upper_half, 1);
//   slap_SetIdentity(lower_half, -1);
// }

enum tiny_ErrorCode tiny_ActiveIneqMask(Matrix* mask, const Matrix dual,
                         const Matrix eval) {
  slap_SetConst(*mask, 0);  // clear before processing
  for (int i = 0; i < eval.rows; ++i) {
    // When variables are on the boundary or violating constraints
    bool active = dual.data[i] > 0 || eval.data[i] > 0;
    slap_SetElement(*mask, i, i, active);
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ProjectOrthantDuals(Matrix* dual, const Matrix new_dual) {
  for (int i = 0; i < dual->rows; ++i) {
    if (new_dual.data[i] > 0) {
      dual->data[i] = new_dual.data[i];
    } else
      dual->data[i] = 0.0;
  }
  return TINY_NO_ERROR;
}

int IsConstrained(tiny_Workspace* work) {
  if (!work->stgs->en_cstr_goal && 
      !work->stgs->en_cstr_inputs && 
      !work->stgs->en_cstr_states) {
    return 0; // unconstrained
  }
  return 1;    
}