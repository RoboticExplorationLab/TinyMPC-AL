#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"

#define NSTATES 2
#define NINPUTS 1
#define NHORIZON 3

class LqrDataTest : public testing::Test {
  public:
    // int nx = 2;
    // int nu = 1;
    // int nh = 2;
    double dt = 0.1;
    double A_data[NSTATES * NSTATES] = {1, 0, 1, 1};  // NOLINT
    double B_data[NSTATES * NINPUTS] = {1, 2};        // NOLINT
    double f_data[NSTATES] = {4, 5};                  // NOLINT
    double t = 1.0;
    double Q_data[NSTATES * NSTATES] = {1, 0, 0, 1};  // NOLINT
    double R_data[NINPUTS * NINPUTS] = {1};           // NOLINT
    double q_data[NSTATES] = {0.1, 0.2};              // NOLINT
    double r_data[NINPUTS] = {-0.6};                  // NOLINT
    double Qf_data[NSTATES * NSTATES] = {0};
    double u_max_data[NINPUTS] = {1.1};
    double u_min_data[NINPUTS] = {-1.1};
    double x_max_data[NSTATES] = {1.6, 1.7};
    double x_min_data[NSTATES] = {-1.6, -1.7};
    double x_ref_data[NSTATES * NHORIZON] = {0.2, 1.1, 2.5, 3.7, 2.1, 4.5};
    double u_ref_data[NINPUTS * (NHORIZON - 1)] = {1, 2};
    double x0_data[NSTATES] = {0.1, 0.2};
    double u0_data[NSTATES] = {1.6};
    double Kd_data[NINPUTS * (NSTATES + 1) * (NHORIZON - 1)] = {0};
    double Pp_data[NSTATES * (NSTATES + 1) * NHORIZON] = {0};
    double reg = 1e-8;
    double input_duals_data[2 * NINPUTS * (NHORIZON - 1)] = {1, 2, 3, 4};
    double state_duals_data[2 * NSTATES * (NHORIZON)] = {1, 2, 3, 4, 5, 6,
                                                        7, 8, 2, 4, 5, 6};
    double goal_duals_data[NSTATES] = {1, 2};
    double reg_min = 1;
    double reg_max = 100;
    double penalty_max = 1e5;
    double penalty_mul = 1;
    int max_primal_iters = 100;
    int max_search_iters = 10;

    const double tol = 1e-8;

  private:
    void SetUp() override {

    }
};

TEST_F(LqrDataTest, InitLtiModelTest) {
  tiny_LtiModel model;
  tiny_InitLtiModel(&model);
  model.nstates = NSTATES;
  model.ninputs = NINPUTS;
  model.dt = dt;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

  EXPECT_EQ(model.dt, dt);
  EXPECT_EQ(model.A.rows, model.nstates);
  EXPECT_EQ(model.A.cols, model.nstates);
  EXPECT_NEAR(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates), 0,
       tol);
  EXPECT_EQ(model.B.rows, model.nstates);
  EXPECT_EQ(model.B.cols, model.ninputs);
  EXPECT_NEAR(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs), 0,
       tol);
  EXPECT_EQ(model.f.rows, model.nstates);
  EXPECT_EQ(model.f.cols, 1);
  EXPECT_NEAR(SumOfSquaredError(f_data, model.f.data, model.nstates), 0, tol);
  EXPECT_EQ(model.x0.rows, model.nstates);
  EXPECT_EQ(model.x0.cols, 1);
  EXPECT_NEAR(SumOfSquaredError(x0_data, model.x0.data, model.nstates), 0, tol);
}

TEST_F(LqrDataTest, KnotPoint) {
  tiny_KnotPoint z;
  tiny_InitKnotPoint(&z);
  tiny_KnotPoint Z[NHORIZON];
  z.dt = dt;
  z.t = t;
  z.x = slap_MatrixFromArray(NSTATES, 1, x0_data);
  z.u = slap_MatrixFromArray(NINPUTS, 1, u0_data);

  EXPECT_EQ(z.dt, dt);
  EXPECT_EQ(z.t, t);
  EXPECT_EQ(z.x.rows, NSTATES);
  EXPECT_EQ(z.x.cols, 1);
  EXPECT_NEAR(SumOfSquaredError(x0_data, z.x.data, z.x.rows), 0, tol);
  EXPECT_EQ(z.u.rows, NINPUTS);
  EXPECT_EQ(z.u.cols, 1);
  EXPECT_NEAR(SumOfSquaredError(u0_data, z.u.data, z.u.rows), 0, tol);

  for (int i = 0; i < NHORIZON; ++i) {
    tiny_InitKnotPoint(&Z[i]);
    Z[i].dt = dt;
    Z[i].t = t + i;
    Z[i].x = slap_MatrixFromArray(NSTATES, 1, x0_data);
    Z[i].u = slap_MatrixFromArray(NINPUTS, 1, u0_data);
  }

  for (int i = 0; i < NHORIZON; ++i) {
    EXPECT_EQ(Z[i].dt, dt);
    EXPECT_EQ(Z[i].t, (t + i));
    EXPECT_EQ(Z[i].x.rows, NSTATES);
    EXPECT_EQ(Z[i].x.cols, 1);
    EXPECT_NEAR(SumOfSquaredError(Z[i].x.data, z.x.data, z.x.rows), 0,  tol);
    EXPECT_EQ(Z[i].u.rows, NINPUTS);
    EXPECT_EQ(Z[i].u.cols, 1);
    EXPECT_NEAR(SumOfSquaredError(Z[i].u.data, z.u.data, z.u.rows), 0, tol);
  }
}

TEST_F(LqrDataTest, Solver) {
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  solver.reg = reg;
  solver.reg_min = reg_min;
  solver.reg_max = reg_max;
  solver.penalty_max = penalty_max;
  solver.penalty_mul = penalty_mul;
  solver.max_primal_iters = max_primal_iters;
  solver.max_search_iters = max_search_iters;

  EXPECT_EQ(solver.reg, reg);
  EXPECT_EQ(solver.reg_max, reg_max);
  EXPECT_EQ(solver.reg_min, reg_min);
  EXPECT_EQ(solver.penalty_max, penalty_max);
  EXPECT_EQ(solver.penalty_mul, penalty_mul);
  EXPECT_EQ(solver.max_primal_iters, max_primal_iters);
  EXPECT_EQ(solver.max_search_iters, max_search_iters);
}

TEST_F(LqrDataTest, ProblemData) {
  const double tol = 1e-8;
  Matrix U_ref[NHORIZON - 1];
  Matrix X_ref[NHORIZON];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix input_duals[NHORIZON - 1];
  Matrix state_duals[NHORIZON];
  double* uptr = u_ref_data;
  double* xptr = x_ref_data;
  double* Kd_ptr = Kd_data;
  double* Pp_ptr = Pp_data;
  double* input_dual_ptr = input_duals_data;
  double* state_dual_ptr = state_duals_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
      uptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd_ptr);
      Kd_ptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, Kd_ptr);
      Kd_ptr += NINPUTS;
      input_duals[i] = slap_MatrixFromArray(NINPUTS, 1, input_dual_ptr);
      input_dual_ptr += NINPUTS;
    }
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pp_ptr);
    Pp_ptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, Pp_ptr);
    Pp_ptr += NSTATES;
    state_duals[i] = slap_MatrixFromArray(NSTATES, 1, state_dual_ptr);
    state_dual_ptr += NSTATES;
  }

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_goal = NSTATES;
  prob.ncstr_inputs = 2 * NINPUTS * (NHORIZON - 1);
  prob.ncstr_states = 2 * NSTATES * NHORIZON;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  prob.q = slap_MatrixFromArray(NSTATES, 1, q_data);
  prob.r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  prob.x_max = slap_MatrixFromArray(NSTATES, 1, x_max_data);
  prob.x_min = slap_MatrixFromArray(NSTATES, 1, x_min_data);
  prob.X_ref = X_ref;
  prob.U_ref = U_ref;
  prob.dt = dt;
  prob.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;
  prob.input_duals = input_duals;
  prob.state_duals = state_duals;
  prob.goal_dual = slap_MatrixFromArray(NSTATES, 1, goal_duals_data);

  EXPECT_EQ(prob.nstates, NSTATES);
  EXPECT_EQ(prob.ninputs, NINPUTS);
  EXPECT_EQ(prob.nhorizon, NHORIZON);
  EXPECT_EQ(prob.ncstr_inputs, 2 * NINPUTS * (NHORIZON - 1));
  EXPECT_EQ(prob.ncstr_states, 2 * NSTATES * NHORIZON);
  EXPECT_EQ(prob.dt, dt);
  EXPECT_NEAR(SumOfSquaredError(prob.Q.data, Q_data, NSTATES * NSTATES), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.R.data, R_data, NINPUTS * NINPUTS), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.q.data, q_data, NSTATES), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.r.data, r_data, NINPUTS), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.Qf.data, Qf_data, NSTATES * NSTATES), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.u_max.data, u_max_data, NINPUTS), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.u_min.data, u_min_data, NINPUTS), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.x_max.data, x_max_data, NSTATES), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.x_min.data, x_min_data, NSTATES), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.x0.data, x0_data, NSTATES), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(prob.goal_dual.data, goal_duals_data, NSTATES), 0, tol);
  uptr = u_ref_data;
  xptr = x_ref_data;
  Kd_ptr = Kd_data;
  Pp_ptr = Pp_data;
  input_dual_ptr = input_duals_data;
  state_dual_ptr = state_duals_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      EXPECT_EQ(prob.U_ref[i].rows, NINPUTS);
      EXPECT_EQ(prob.U_ref[i].cols, 1);
      EXPECT_NEAR(SumOfSquaredError(prob.U_ref[i].data, uptr, NINPUTS), 0, tol);
      uptr += NINPUTS;
      EXPECT_EQ(prob.K[i].rows, NINPUTS);
      EXPECT_EQ(prob.K[i].cols, NSTATES);
      EXPECT_NEAR(SumOfSquaredError(prob.K[i].data, Kd_ptr, NINPUTS * NSTATES), 0, tol);
      Kd_ptr += NINPUTS * NSTATES;
      EXPECT_EQ(prob.d[i].rows, NINPUTS);
      EXPECT_EQ(prob.d[i].cols, 1);
      EXPECT_NEAR(SumOfSquaredError(prob.d[i].data, Kd_ptr, NINPUTS), 0, tol);
      Kd_ptr += NINPUTS;
      EXPECT_EQ(prob.input_duals[i].rows, NINPUTS);
      EXPECT_EQ(prob.input_duals[i].cols, 1);
      EXPECT_NEAR(SumOfSquaredError(prob.input_duals[i].data, input_dual_ptr,
                             NINPUTS), 0, tol);
      input_dual_ptr += NINPUTS;
    }
    EXPECT_EQ(prob.X_ref[i].rows, NSTATES);
    EXPECT_EQ(prob.X_ref[i].cols, 1);
    EXPECT_NEAR(SumOfSquaredError(prob.X_ref[i].data, xptr, NSTATES), 0, tol);
    xptr += NSTATES;
    EXPECT_EQ(prob.P[i].rows, NSTATES);
    EXPECT_EQ(prob.P[i].cols, NSTATES);
    EXPECT_NEAR(SumOfSquaredError(prob.P[i].data, Pp_ptr, NSTATES * NSTATES), 0, tol);
    Pp_ptr += NSTATES * NSTATES;
    EXPECT_EQ(prob.p[i].rows, NSTATES);
    EXPECT_EQ(prob.p[i].cols, 1);
    EXPECT_NEAR(SumOfSquaredError(prob.p[i].data, Pp_ptr, NSTATES), 0, tol);
    Pp_ptr += NSTATES;
    EXPECT_EQ(prob.state_duals[i].rows, NSTATES);
    EXPECT_EQ(prob.state_duals[i].cols, 1);
    EXPECT_NEAR(SumOfSquaredError(prob.state_duals[i].data, state_dual_ptr, NSTATES), 0,
         tol);
    state_dual_ptr += NSTATES;
  }
}