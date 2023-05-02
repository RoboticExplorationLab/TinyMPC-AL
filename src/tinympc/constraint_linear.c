#include "constraint_linear.h"

sfloat tiny_RiccatiConvergence(const tiny_ProblemData prob) {
  sfloat norm_d_max = 0.0;
  for (int k = 0; k < prob.nhorizon - 1; ++k) {
    sfloat norm_d = slap_NormTwo(prob.d[k]);
    if (norm_d > norm_d_max) {
      norm_d_max = norm_d;
    }
  }
  return norm_d_max;
}

// // [u-p.u_max;-u + p.u_min]
// void tiny_IneqInputs(Matrix* cu, const tiny_ProblemData prob,
//                      const Matrix u) {
//   slap_SetConst(*cu, 0);  // clear before processing
//   Matrix upper_half = slap_CreateSubMatrix(*cu, 0, 0, prob.ninputs,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cu, prob.ninputs, 0, prob.ninputs, 1);
//   slap_MatrixAddition(upper_half, u, prob.u_max, -1);
//   slap_MatrixAddition(lower_half, prob.u_min, u, -1);
// }

// Input constraints: A*u - b <= 0
void tiny_IneqInputs(Matrix* cu, const tiny_ProblemData prob,
                     const Matrix u) {
  // slap_SetConst(*cu, 0);  // clear before processing
  slap_MatMulAdd(*cu, prob.Acu, u, 1.0, 0.0);
  slap_MatrixAddition(*cu, *cu, prob.bcu, -1.0);
}

// // [u_max, -u_min]
// void tiny_IneqInputsOffset(Matrix* cu, const tiny_ProblemData prob) {
//   Matrix upper_half = slap_CreateSubMatrix(*cu, 0, 0, prob.ninputs,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cu, prob.ninputs, 0, prob.ninputs, 1);
//   slap_Copy(upper_half, prob.u_max);
//   slap_Copy(lower_half, prob.u_min);
//   slap_ScaleByConst(lower_half, -1);
// }

// void tiny_IneqInputsJacobian(Matrix* ineq_jac, const tiny_ProblemData prob) {
//   slap_SetConst(*ineq_jac, 0);  // clear before processing
//   Matrix upper_half =
//       slap_CreateSubMatrix(*ineq_jac, 0, 0, prob.ninputs, prob.ninputs);
//   Matrix lower_half = slap_CreateSubMatrix(*ineq_jac, prob.ninputs, 0,
//                                            prob.ninputs, prob.ninputs);
//   slap_SetIdentity(upper_half, 1);
//   slap_SetIdentity(lower_half, -1);
// }

// // [x-p.x_max;-x + p.x_min]
// void tiny_IneqStates(Matrix* cx, const tiny_ProblemData prob,
//                      const Matrix x) {
//   slap_SetConst(*cx, 0);  // clear before processing
//   Matrix upper_half = slap_CreateSubMatrix(*cx, 0, 0, prob.nstates,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cx, prob.nstates, 0, prob.nstates, 1);
//   slap_MatrixAddition(upper_half, x, prob.x_max, -1);
//   slap_MatrixAddition(lower_half, prob.x_min, x, -1);
// }

// State constraints: A*x - b <= 0
void tiny_IneqStates(Matrix* cx, const tiny_ProblemData prob,
                     const Matrix x) {
  // slap_SetConst(*cu, 0);  // clear before processing
  slap_MatMulAdd(*cx, prob.Acx, x, 1.0, 0.0);
  slap_MatrixAddition(*cx, *cx, prob.bcx, -1.0);
}

// // [x_max, -x_min]
// void tiny_IneqStatesOffset(Matrix* cx, const tiny_ProblemData prob) {
//   Matrix upper_half = slap_CreateSubMatrix(*cx, 0, 0, prob.nstates,
//   1); Matrix lower_half =
//       slap_CreateSubMatrix(*cx, prob.nstates, 0, prob.nstates, 1);
//   slap_Copy(upper_half, prob.x_max);
//   slap_Copy(lower_half, prob.x_min);
//   slap_ScaleByConst(lower_half, -1);
// }

// void tiny_IneqStatesJacobian(Matrix* ineq_jac, const tiny_ProblemData prob) {
//   slap_SetConst(*ineq_jac, 0);  // clear before processing
//   Matrix upper_half =
//       slap_CreateSubMatrix(*ineq_jac, 0, 0, prob.nstates, prob.nstates);
//   Matrix lower_half = slap_CreateSubMatrix(*ineq_jac, prob.nstates, 0,
//                                            prob.nstates, prob.nstates);
//   slap_SetIdentity(upper_half, 1);
//   slap_SetIdentity(lower_half, -1);
// }

void tiny_ActiveIneqMask(Matrix* mask, const Matrix dual, const Matrix ineq) {
  slap_SetConst(*mask, 0);  // clear before processing
  for (int i = 0; i < ineq.rows; ++i) {
    // When variables are on the boundary or violating constraints
    bool active = dual.data[i] > 0 || ineq.data[i] > 0;
    slap_SetElement(*mask, i, i, active);
  }
}

void tiny_ClampIneqDuals(Matrix* dual, const Matrix new_dual) {
  for (int i = 0; i < dual->rows; ++i) {
    if (new_dual.data[i] > 0) {
      dual->data[i] = new_dual.data[i];
    } else
      dual->data[i] = 0.0;
  }
}