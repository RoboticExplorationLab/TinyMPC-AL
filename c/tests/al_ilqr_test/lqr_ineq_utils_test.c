// Test all inequality-related functions and tiny_RiccatiConvergence

#include "constrained_ilqr.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

void IneqInputsTest() {
  const sfloat tol = 1e-6;
  sfloat u_max_data[NINPUTS] = {1, 1};
  sfloat u_min_data[NINPUTS] = {-1, -1};
  sfloat u_data[NINPUTS] = {1.1, 0.8};
  sfloat ans[NINPUTS * 2] = {
      u_data[0] - u_max_data[0], u_data[1] - u_max_data[1],
      -u_data[0] + u_min_data[0], -u_data[1] + u_min_data[1]};
  sfloat ineq_data[NINPUTS * 2];

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2 * NSTATES;
  prob.ncstr_inputs = 2 * NINPUTS;
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);
  Matrix mat = slap_MatrixFromArray(NINPUTS * 2, 1, ineq_data);
  tiny_IneqInputs(&mat, prob, u);
  // slap_PrintMatrix(mat);
  TEST(SumOfSquaredError(mat.data, ans, NINPUTS * 2) < tol);
}

void IneqInputsJacobianTest() {
  const sfloat tol = 1e-6;
  sfloat u_max_data[NINPUTS] = {1, 1};
  sfloat u_min_data[NINPUTS] = {-1, -1};
  sfloat ans[NINPUTS * 2 * NINPUTS] = {1, 0, -1, 0, 0, 1, 0, -1};
  sfloat jac_data[NINPUTS * 2 * NINPUTS];

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2 * NSTATES;
  prob.ncstr_inputs = 2 * NINPUTS;
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix mat = slap_MatrixFromArray(NINPUTS * 2, NINPUTS, jac_data);
  tiny_IneqInputsJacobian(&mat, prob);
  // slap_PrintMatrix(mat);
  TEST(SumOfSquaredError(mat.data, ans, NINPUTS * 2 * NINPUTS) < tol);
}

void ActiveIneqMaskTest() {
  const sfloat tol = 1e-6;
  sfloat u_max_data[NINPUTS] = {1, 1};
  sfloat u_min_data[NINPUTS] = {-1, -1};
  sfloat u_data[NINPUTS] = {-2, 2};
  sfloat ans1[NINPUTS * 2] = {u_data[0] - u_max_data[0],    //-3
                              u_data[1] - u_max_data[1],    // 1
                              -u_data[0] + u_min_data[0],   // 1
                              -u_data[1] + u_min_data[1]};  // -3
  sfloat dual_data[NINPUTS * 2] = {1, 0, 2, 0};
  sfloat ans2[NINPUTS * 2 * NINPUTS * 2] = {1, 0, 0, 0, 0, 1, 0, 0,
                                            0, 0, 1, 0, 0, 0, 0, 0};
  sfloat ineq_data[NINPUTS * 2];
  sfloat mask_data[NINPUTS * 2 * NINPUTS * 2];

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2 * NSTATES;
  prob.ncstr_inputs = 2 * NINPUTS;
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);
  Matrix ineq = slap_MatrixFromArray(NINPUTS * 2, 1, ineq_data);
  tiny_IneqInputs(&ineq, prob, u);
  Matrix input_dual = slap_MatrixFromArray(NINPUTS * 2, 1, dual_data);
  Matrix mask = slap_MatrixFromArray(NINPUTS * 2, NINPUTS * 2, mask_data);
  tiny_ActiveIneqMask(&mask, input_dual, ineq);
  // slap_PrintMatrix(mask);
  TEST(SumOfSquaredError(ineq.data, ans1, NINPUTS * 2) < tol);
  TEST(SumOfSquaredError(mask.data, ans2, NINPUTS * 2 * NINPUTS * 2) < tol);
}
void RiccatiConvergenceTest() {
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nhorizon = 3;
  prob.ninputs = 2;
  sfloat d_data[4] = {1.2, -0.3, -2.1, 3.1};
  sfloat ans = 3.744329045369811;
  Matrix d[2];
  sfloat* dptr = d_data;
  for (int k = 0; k < prob.nhorizon - 1; ++k) {
    d[k] = slap_MatrixFromArray(prob.ninputs, 1, dptr);
    dptr += prob.ninputs;
  }
  prob.d = d;
  sfloat norm_d_max = tiny_RiccatiConvergence(prob);
  TESTAPPROX(norm_d_max, ans, 1e-6);
}

void ClampIneqDualsTest() {
  sfloat dual_data[2] = {-2, -1};
  sfloat new_dual_data[2] = {10, -10};
  sfloat ans[2] = {10, 0};
  Matrix dual = slap_MatrixFromArray(2, 1, dual_data);
  Matrix new_dual = slap_MatrixFromArray(2, 1, new_dual_data);
  tiny_ClampIneqDuals(&dual, new_dual);
  TEST(SumOfSquaredError(dual.data, ans, 2) < 1e-8);
}

int main() {
  IneqInputsTest();
  IneqInputsJacobianTest();
  ActiveIneqMaskTest();
  RiccatiConvergenceTest();
  ClampIneqDualsTest();
  PrintTestResult();
  return TestResult();
}
