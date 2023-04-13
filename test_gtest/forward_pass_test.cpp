#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"
#include "data/forward_pass_data.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

class ForwardPassTest : public testing::Test {
  public:
    double A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                        0.1, 0, 1, 0, 0, 0.1, 0, 1};
    double B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
    double f_data[NSTATES] = {0, 0, 0, 0};
    // double x0_data[NSTATES] = {5,7,2,-1.4};

    const double tol = 1e-8;

    Matrix X[NHORIZON];
    Matrix Xsln[NHORIZON];
    Matrix U[NHORIZON - 1];
    Matrix K[NHORIZON - 1];
    Matrix d[NHORIZON - 1];

    tiny_LtiModel model;
    tiny_ProblemData prob;

    double* xptr = x_data;
    double* xsol_ptr = xsol_data;
    double* uptr = u_data;
    double* Kptr = K_data;
    double* dptr = d_data;

  private:
    void SetUp() override {
      tiny_InitLtiModel(&model);
      tiny_InitProblemData(&prob);

      model.dt = 0.1;
      model.ninputs = NSTATES;
      model.nstates = NINPUTS;
      model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
      model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
      model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
      // model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

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
    }
};


TEST_F(ForwardPassTest, SetUpTest) {
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      EXPECT_EQ(U[i].rows, NINPUTS);
      EXPECT_EQ(U[i].cols, 1);
      EXPECT_NEAR(SumOfSquaredError(U[i].data, uptr, NINPUTS), 0, tol);
      uptr += NINPUTS;
      EXPECT_EQ(prob.K[i].rows, NINPUTS);
      EXPECT_EQ(prob.K[i].cols, NSTATES);
      EXPECT_NEAR(SumOfSquaredError(prob.K[i].data, Kptr, NINPUTS * NSTATES), 0, tol);
      Kptr += NINPUTS * NSTATES;
      EXPECT_EQ(prob.d[i].rows, NINPUTS);
      EXPECT_EQ(prob.d[i].cols, 1);
      EXPECT_NEAR(SumOfSquaredError(prob.d[i].data, dptr, NINPUTS), 0, tol);
      dptr += NINPUTS;
    }
    EXPECT_EQ(X[i].rows, NSTATES);
    EXPECT_EQ(X[i].cols, 1);
    EXPECT_NEAR(SumOfSquaredError(X[i].data, xptr, NSTATES), 0, tol);
    xptr += NSTATES;
    EXPECT_EQ(Xsln[i].rows, NSTATES);
    EXPECT_EQ(Xsln[i].cols, 1);
    EXPECT_NEAR(SumOfSquaredError(Xsln[i].data, xsol_ptr, NSTATES), 0, tol);
    xsol_ptr += NSTATES;
  }
}

TEST_F(ForwardPassTest, DiscreteDynamics) {
  // Include discrete dynamics test
  tiny_ForwardPassLti(X, U, prob, model);
  for (int i = 0; i < NHORIZON; ++i) {
    EXPECT_NEAR(SumOfSquaredError(X[i].data, Xsln[i].data, NSTATES), 0, tol);
  }
}