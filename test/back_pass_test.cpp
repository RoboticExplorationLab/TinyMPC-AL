#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"
#include "data/back_pass_data.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

class BackPassTest : public testing::Test {
  public:
    double A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                        0.1, 0, 1, 0, 0, 0.1, 0, 1};
    double B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
    double f_data[NSTATES] = {0};
    double x0_data[NSTATES] = {5, 7, 2, -1.4};
    double Xref_data[NSTATES * NHORIZON] = {0};
    double Uref_data[NINPUTS * (NHORIZON - 1)] = {0};
    // double X_data[NSTATES*NHORIZON] = {0};
    // double U_data[NINPUTS*(NHORIZON-1)] = {0};
    double K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
    double d_data[NINPUTS * (NHORIZON - 1)] = {0};
    double P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
    double p_data[NSTATES * NHORIZON] = {0};
    double Q_data[NSTATES * NSTATES] = {0};
    double R_data[NINPUTS * NINPUTS] = {0};
    double Qf_data[NSTATES * NSTATES] = {0};
    double umin_data[NINPUTS] = {-2, -2};
    double umax_data[NINPUTS] = {2, 2};
    const double tol = 1e-8;

    tiny_LtiModel model;
    tiny_ProblemData prob;
    tiny_Solver solver;

    Matrix Xref[NHORIZON];
    Matrix Uref[NHORIZON - 1];
    Matrix K[NHORIZON - 1];
    Matrix d[NHORIZON - 1];
    Matrix P[NHORIZON];
    Matrix p[NHORIZON];
    
  private:
    void SetUp() override {
      tiny_InitLtiModel(&model);
      tiny_InitProblemData(&prob);
      tiny_InitSolver(&solver);

      model.ninputs = NSTATES;
      model.nstates = NINPUTS;
      model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
      model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
      model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
      model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
    }
};

// U, X, Psln
TEST_F(BackPassTest, BackPass) {
  double* Xref_ptr = Xref_data;
  double* Uref_ptr = Uref_data;
  double* Kptr = K_data;
  double* dptr = d_data;
  double* Pptr = P_data;
  double* pptr = p_data;
  
  
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    Xref_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
  }

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 2 * NINPUTS * (NHORIZON - 1);
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 100 * 1e-1);
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, umax_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, umin_data);
  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  solver.reg = 1e-8;
  solver.penalty_mul = 10;

  double G_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix G_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       G_temp_data);
  tiny_BackwardPassLti(&prob, solver, model, &G_temp);
  EXPECT_NEAR(SumOfSquaredError(d_data, dsln_data, (NHORIZON - 1) * NINPUTS), tol);
}