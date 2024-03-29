// Task: Test AL-LQR on double integrator with input/state box constraints and
// goal constraint. Scenerio: drive from initial state to goal state.

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51
// U, X, Psln

class MpcLtiTest : public testing::Test {
 public:
  double A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                      0.1, 0, 1, 0, 0, 0.1, 0, 1};
  double B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
  double f_data[NSTATES] = {0};
  double x0_data[NSTATES] = {5, 7, 2, -1.4};
  double xg_data[NSTATES] = {0};
  double Xref_data[NSTATES * NHORIZON] = {0};
  double Uref_data[NINPUTS * (NHORIZON - 1)] = {0};
  double X_data[NSTATES * NHORIZON] = {0};
  double U_data[NINPUTS * (NHORIZON - 1)] = {0};
  double K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  double d_data[NINPUTS * (NHORIZON - 1)] = {0};
  double P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  double p_data[NSTATES * NHORIZON] = {0};
  double Q_data[NSTATES * NSTATES] = {0};
  double R_data[NINPUTS * NINPUTS] = {0};
  double Qf_data[NSTATES * NSTATES] = {0};
  double umin_data[NINPUTS] = {-6, -6};
  double umax_data[NINPUTS] = {6, 6};
  double xmin_data[NSTATES] = {-2, -2, -2, -2};
  double xmax_data[NSTATES] = {6, 8, 3, 2};
  double input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0};
  double state_dual_data[2 * NSTATES * (NHORIZON)] = {0};
  double goal_dual_data[NSTATES] = {0};

  tiny_LtiModel model;
  tiny_ProblemData prob;
  tiny_Solver solver;

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix input_duals[NHORIZON - 1];
  Matrix state_duals[NHORIZON];

  Matrix xg;

 protected:
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
    xg = slap_MatrixFromArray(NSTATES, 1, xg_data);
  }
};

TEST_F(MpcLtiTest, DoubleIntegrator) {
  double* Xptr = X_data;
  double* Xref_ptr = Xref_data;
  double* Uptr = U_data;
  double* Uref_ptr = Uref_data;
  double* Kptr = K_data;
  double* dptr = d_data;
  double* Pptr = P_data;
  double* pptr = p_data;
  double* udual_ptr = input_dual_data;
  double* xdual_ptr = state_dual_data;
  
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_MatrixFromArray(NINPUTS, 1, Uptr);
      slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
      input_duals[i] = slap_MatrixFromArray(2 * NINPUTS, 1, udual_ptr);
      udual_ptr += 2 * NINPUTS;
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    slap_Copy(Xref[i], xg);
    Xref_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
    state_duals[i] = slap_MatrixFromArray(2 * NSTATES, 1, xdual_ptr);
    xdual_ptr += 2 * NSTATES;
  }
  
  slap_Copy(X[0], model.x0);
  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 2 * NINPUTS;
  prob.ncstr_states = 2 * NSTATES;
  prob.ncstr_goal = NSTATES;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 100 * 1e-1);
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, umax_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, umin_data);
  prob.x_max = slap_MatrixFromArray(NSTATES, 1, xmax_data);
  prob.x_min = slap_MatrixFromArray(NSTATES, 1, xmin_data);
  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;
  prob.input_duals = input_duals;
  prob.state_duals = state_duals;
  prob.goal_dual = slap_MatrixFromArray(NSTATES, 1, goal_dual_data);

  solver.max_primal_iters = 10;
  tiny_MpcLti(X, U, &prob, &solver, model, 0);

  for (int k = 0; k < NHORIZON - 1; ++k) {
    // tiny_Print(X[k]);
    EXPECT_LT(slap_NormInf(U[k]), slap_NormInf(prob.u_max) + solver.cstr_tol);
    for (int i = 0; i < NSTATES; ++i) {
      EXPECT_LT(X[k].data[i], xmax_data[i] + solver.cstr_tol);
      EXPECT_GT(X[k].data[i], xmin_data[i] - solver.cstr_tol);
    }
  }
  
  EXPECT_LT(SumOfSquaredError(X[NHORIZON - 1].data, xg_data, NSTATES), solver.cstr_tol);
}
