#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "augmented_lagrangian_lqr.h"
#include "simpletest/simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

void IneqInputsTest() {
  const double tol = 1e-8;
  double u_max_data[NINPUTS] = {1, 1};
  double u_min_data[NINPUTS] = {-1, -1};
  double u_data[NINPUTS] = {1.1, 0.8};
  double ans[NINPUTS*2] = {u_data[0] - u_max_data[0],
                           u_data[1] - u_max_data[1],
                           -u_data[0] + u_min_data[0],
                           -u_data[1] + u_min_data[1]};
  double ineq_data[NINPUTS*2];

  tiny_ProblemData prob = kDefaultProblemData;
  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2*NSTATES;
  prob.ncstr_inputs = 2*NINPUTS;
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);
  Matrix mat = slap_MatrixFromArray(NINPUTS*2, 1, ineq_data);
  tiny_IneqInputs(mat, prob, u);
  // slap_PrintMatrix(mat);
  TEST(SumOfSquaredError(mat.data, ans, NINPUTS*2) < tol);
}

void IneqInputsJacobianTest() {
  const double tol = 1e-8;
  double u_max_data[NINPUTS] = {1, 1};
  double u_min_data[NINPUTS] = {-1, -1};
  double u_data[NINPUTS] = {1.1, 0.8};
  double ans[NINPUTS*2*NINPUTS] = {1, 0, -1, 0, 0, 1, 0, -1};
  double jac_data[NINPUTS*2*NINPUTS];

  tiny_ProblemData prob = kDefaultProblemData;
  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2*NSTATES;
  prob.ncstr_inputs = 2*NINPUTS;
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);
  Matrix mat = slap_MatrixFromArray(NINPUTS*2, NINPUTS, jac_data);
  tiny_IneqInputsJacobian(mat, prob, u);
  // slap_PrintMatrix(mat);
  TEST(SumOfSquaredError(mat.data, ans, NINPUTS*2*NINPUTS) < tol);
}

void ActiveIneqMaskTest() {
  const double tol = 1e-8;
  double u_max_data[NINPUTS] = {1, 1};
  double u_min_data[NINPUTS] = {-1, -1};
  double u_data[NINPUTS] = {-2, 2};
  double ans1[NINPUTS*2] = {u_data[0] - u_max_data[0], //-3
                           u_data[1] - u_max_data[1], // 1
                           -u_data[0] + u_min_data[0],  // 1
                           -u_data[1] + u_min_data[1]}; // -3
  double dual_data[NINPUTS*2] = {1, 0, 2, 0};
  double ans2[NINPUTS*2*NINPUTS*2] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,0};
  double ineq_data[NINPUTS*2];
  double mask_data[NINPUTS*2*NINPUTS*2];

  tiny_ProblemData prob = kDefaultProblemData;
  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_states = 2*NSTATES;
  prob.ncstr_inputs = 2*NINPUTS;
  prob.ncstr_goal = NSTATES;
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);
  Matrix ineq = slap_MatrixFromArray(NINPUTS*2, 1, ineq_data);
  tiny_IneqInputs(ineq, prob, u);
  Matrix input_dual = slap_MatrixFromArray(NINPUTS*2, 1, dual_data);
  Matrix mask = slap_MatrixFromArray(NINPUTS*2, NINPUTS*2, mask_data);
  tiny_ActiveIneqMask(mask, input_dual, ineq);
  // slap_PrintMatrix(mask);
  TEST(SumOfSquaredError(ineq.data, ans1, NINPUTS*2) < tol);
  TEST(SumOfSquaredError(mask.data, ans2, NINPUTS*2*NINPUTS*2) < tol);
}

int main() {
  IneqInputsTest();
  IneqInputsJacobianTest();
  ActiveIneqMaskTest();
  PrintTestResult();
  return TestResult();
}
