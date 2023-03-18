#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data/forward_pass_data.h"
#include "gtest/gtest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/dynamics_lti.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

double A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                    0.1, 0, 1, 0, 0, 0.1, 0, 1};
double B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
double f_data[NSTATES] = {0, 0, 0, 0};
// double x0_data[NSTATES] = {5,7,2,-1.4};

void ForwardPassTest() {
  const double tol = 1e-8;
  tiny_LtiModel model;
  tiny_InitLtiModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  model.dt = 0.1;
  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
  // model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

  Matrix X[NHORIZON];
  Matrix Xsln[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];

  double* xptr = x_data;
  double* xsol_ptr = xsol_data;
  double* uptr = u_data;
  double* Kptr = K_data;
  double* dptr = d_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
      uptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
    Xsln[i] = slap_MatrixFromArray(NSTATES, 1, xsol_ptr);
    xsol_ptr += NSTATES;
  }

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.x0 = X[0];  // check if possible
  prob.K = K;
  prob.d = d;

  uptr = u_data;
  xsol_ptr = xsol_data;
  xptr = x_data;
  Kptr = K_data;
  dptr = d_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      TEST(U[i].rows == NINPUTS);
      TEST(U[i].cols == 1);
      TEST(SumOfSquaredError(U[i].data, uptr, NINPUTS) < tol);
      uptr += NINPUTS;
      TEST(prob.K[i].rows == NINPUTS);
      TEST(prob.K[i].cols == NSTATES);
      TEST(SumOfSquaredError(prob.K[i].data, Kptr, NINPUTS * NSTATES) < tol);
      Kptr += NINPUTS * NSTATES;
      TEST(prob.d[i].rows == NINPUTS);
      TEST(prob.d[i].cols == 1);
      TEST(SumOfSquaredError(prob.d[i].data, dptr, NINPUTS) < tol);
      dptr += NINPUTS;
    }
    TEST(X[i].rows == NSTATES);
    TEST(X[i].cols == 1);
    TEST(SumOfSquaredError(X[i].data, xptr, NSTATES) < tol);
    xptr += NSTATES;
    TEST(Xsln[i].rows == NSTATES);
    TEST(Xsln[i].cols == 1);
    TEST(SumOfSquaredError(Xsln[i].data, xsol_ptr, NSTATES) < tol);
    xsol_ptr += NSTATES;
  }

  // Include discrete dynamics test
  tiny_ForwardPassLti(X, U, prob, model);
  for (int i = 0; i < NHORIZON; ++i) {
    TEST(SumOfSquaredError(X[i].data, Xsln[i].data, NSTATES) < tol);
  }
}

int main() {
  ForwardPassTest();
  PrintTestResult();
  return TestResult();
}
