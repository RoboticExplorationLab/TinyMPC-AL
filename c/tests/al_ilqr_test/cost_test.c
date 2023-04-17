#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constrained_ilqr.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 2

sfloat x_data[NSTATES] = {1.1, 1.2, 1.3, -4.3};
sfloat u_data[NSTATES] = {-2.1, 1.1};
sfloat x_ref_data[NSTATES * NHORIZON] = {1.1, 1.2, 1.3, -4.2,
                                         1.2, 1.3, 1.3, -4.3};
sfloat u_ref_data[NINPUTS * NHORIZON] = {-2.1, 1.4, -2.2, 1.5};
sfloat Q_data[NSTATES * NSTATES] = {0};  // NOLINT
sfloat R_data[NINPUTS * NINPUTS] = {0};  // NOLINT
sfloat q_data[NSTATES] = {0};            // NOLINT
sfloat r_data[NINPUTS] = {0};            // NOLINT
sfloat Qf_data[NSTATES * NSTATES] = {0};
sfloat ans_stage[2] = {0.04549999999999994, 0.1314999999999999};
sfloat ans_term = 0.0049999999999999975;
sfloat ans_gradx[NSTATES] = {0.0, 0.0, 0.0, -0.009999999999999966};
sfloat ans_gradu[NINPUTS] = {0.0, -0.2999999999999998};
sfloat ans_gradxf[NSTATES] = {-0.04999999999999993, -0.050000000000000044, 0.0,
                              0.0};

void AddCostTest() {
  const sfloat tol = 1e-6;
  sfloat cost = 0;
  Matrix U_ref[NHORIZON];
  Matrix X_ref[NHORIZON];
  sfloat* uptr = u_ref_data;
  sfloat* xptr = x_ref_data;
  for (int i = 0; i < NHORIZON; ++i) {
    U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
    uptr += NINPUTS;
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
  }
  Matrix x = slap_MatrixFromArray(NSTATES, 1, x_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 0.1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1);
  prob.q = slap_MatrixFromArray(NSTATES, 1, q_data);
  slap_SetConst(prob.q, 1);
  prob.r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  slap_SetConst(prob.r, 2);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 0.5);
  prob.X_ref = X_ref;
  prob.U_ref = U_ref;

  for (int k = 0; k < 2; ++k) {
    tiny_AddStageCost(&cost, prob, x, u, k);
    TESTAPPROX(cost, ans_stage[k], tol);
  }
  cost = 0;
  tiny_AddTerminalCost(&cost, prob, x);
  TESTAPPROX(cost, ans_term, tol);
}

void ExpandCostTest() {
  const sfloat tol = 1e-6;
  Matrix U_ref[NHORIZON];
  Matrix X_ref[NHORIZON];
  sfloat* uptr = u_ref_data;
  sfloat* xptr = x_ref_data;
  for (int i = 0; i < NHORIZON; ++i) {
    U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
    uptr += NINPUTS;
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
  }
  Matrix x = slap_MatrixFromArray(NSTATES, 1, x_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 0.1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1);
  prob.q = slap_MatrixFromArray(NSTATES, 1, q_data);
  slap_SetConst(prob.q, 1);
  prob.r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  slap_SetConst(prob.r, 2);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 0.5);
  prob.X_ref = X_ref;
  prob.U_ref = U_ref;

  sfloat hessx_data[NSTATES * NSTATES];
  sfloat gradx_data[NSTATES];
  sfloat hessu_data[NSTATES * NSTATES];
  sfloat gradu_data[NSTATES];
  Matrix hessx = slap_MatrixFromArray(NSTATES, NSTATES, hessx_data);
  Matrix gradx = slap_MatrixFromArray(NSTATES, 1, gradx_data);
  Matrix hessu = slap_MatrixFromArray(NINPUTS, NINPUTS, hessu_data);
  Matrix gradu = slap_MatrixFromArray(NINPUTS, 1, gradu_data);
  tiny_ExpandStageCost(&hessx, &gradx, &hessu, &gradu, prob, x, u, 0);
  TEST(SumOfSquaredError(hessx.data, Q_data, NSTATES * NSTATES) < tol);
  TEST(SumOfSquaredError(hessu.data, R_data, NSTATES) < tol);
  TEST(SumOfSquaredError(gradx.data, ans_gradx, NINPUTS * NINPUTS) < tol);
  TEST(SumOfSquaredError(gradu.data, ans_gradu, NINPUTS) < tol);
  tiny_ExpandTerminalCost(&hessx, &gradx, prob, x);
  TEST(SumOfSquaredError(hessx.data, Qf_data, NSTATES * NSTATES) < tol);
  TEST(SumOfSquaredError(gradx.data, ans_gradxf, NSTATES) < tol);
}

int main() {
  AddCostTest();
  ExpandCostTest();
  PrintTestResult();
  return TestResult();
}
