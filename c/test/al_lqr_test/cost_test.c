#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constrained_lqr.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 2

double x_data[NSTATES] = {1.1, 1.2, 1.3, -4.3};
double u_data[NSTATES] = {-2.1, 1.1};
double x_ref_data[NSTATES * NHORIZON] = {1.1, 1.2, 1.3, -4.2,
                                         1.2, 1.3, 1.3, -4.3};
double u_ref_data[NINPUTS * NHORIZON] = {-2.1, 1.4, -2.2, 1.5};
double Q_data[NSTATES * NSTATES] = {0};  // NOLINT
double R_data[NINPUTS * NINPUTS] = {0};  // NOLINT
double q_data[NSTATES] = {0};            // NOLINT
double r_data[NINPUTS] = {0};            // NOLINT
double Qf_data[NSTATES * NSTATES] = {0};
double ans_stage[2] = {0.04549999999999994, 0.1314999999999999};
double ans_term = 0.0049999999999999975;
double ans_gradx[NSTATES] = {-0.11, -0.12, -0.13,  0.42};
double ans_gradu[NINPUTS] = {0.21, -0.14};
double ans_gradxf[NSTATES] = {-0.6, -0.65, -0.65,  2.15};

void AddCostTest() {
  const double tol = 1e-8;
  double cost = 0;
  Matrix U_ref[NHORIZON];
  Matrix X_ref[NHORIZON];
  double* uptr = u_ref_data;
  double* xptr = x_ref_data;
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
  const double tol = 1e-8;
  Matrix U_ref[NHORIZON];
  Matrix X_ref[NHORIZON];
  double* uptr = u_ref_data;
  double* xptr = x_ref_data;
  for (int i = 0; i < NHORIZON; ++i) {
    U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
    uptr += NINPUTS;
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
  }

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 0.1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 0.1);
  prob.q = slap_MatrixFromArray(NSTATES, 1, q_data);
  slap_SetConst(prob.q, 1);
  prob.r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  slap_SetConst(prob.r, 2);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 0.5);
  prob.X_ref = X_ref;
  prob.U_ref = U_ref;

  double hessx_data[NSTATES * NSTATES];
  double gradx_data[NSTATES];
  double hessu_data[NSTATES * NSTATES];
  double gradu_data[NSTATES];
  Matrix hessx = slap_MatrixFromArray(NSTATES, NSTATES, hessx_data);
  Matrix gradx = slap_MatrixFromArray(NSTATES, 1, gradx_data);
  Matrix hessu = slap_MatrixFromArray(NINPUTS, NINPUTS, hessu_data);
  Matrix gradu = slap_MatrixFromArray(NINPUTS, 1, gradu_data);
  tiny_ExpandStageCost(&hessx, &gradx, &hessu, &gradu, prob, 0);
  TEST(SumOfSquaredError(hessx.data, Q_data, NSTATES * NSTATES) < tol);
  TEST(SumOfSquaredError(hessu.data, R_data, NSTATES) < tol);
  TEST(SumOfSquaredError(gradx.data, ans_gradx, NINPUTS * NINPUTS) < tol);
  TEST(SumOfSquaredError(gradu.data, ans_gradu, NINPUTS) < tol);
  tiny_ExpandTerminalCost(&hessx, &gradx, prob);
  TEST(SumOfSquaredError(hessx.data, Qf_data, NSTATES * NSTATES) < tol);
  TEST(SumOfSquaredError(gradx.data, ans_gradxf, NSTATES) < tol);
}

int main() {
  AddCostTest();
  ExpandCostTest();
  PrintTestResult();
  return TestResult();
}