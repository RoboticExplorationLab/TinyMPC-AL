// Task: Test AL-LQR on double integrator with input/state box constraints and
// goal constraint. Scenerio: drive from initial state to goal state.

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/mpc_lti.h"
#include "tinympc/utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51
// U, X, Psln
void MpcLtiTest() {
  // sfloat tol = 1e-4;
  sfloat A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                      0.1, 0, 1, 0, 0, 0.1, 0, 1};
  sfloat B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
  sfloat f_data[NSTATES] = {0};
  sfloat x0_data[NSTATES] = {5, 7, 2, -1.4};
  sfloat xg_data[NSTATES] = {0};
  sfloat Xref_data[NSTATES * NHORIZON] = {0};
  sfloat Uref_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat X_data[NSTATES * NHORIZON] = {0};
  sfloat U_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  sfloat d_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  sfloat p_data[NSTATES * NHORIZON] = {0};
  sfloat Q_data[NSTATES * NSTATES] = {0};
  sfloat R_data[NINPUTS * NINPUTS] = {0};
  sfloat Qf_data[NSTATES * NSTATES] = {0};
  sfloat umin_data[NINPUTS] = {-5, -5};
  sfloat umax_data[NINPUTS] = {5, 5};
  sfloat xmin_data[NSTATES] = {-2, -2, -2, -2};
  sfloat xmax_data[NSTATES] = {6, 8, 3, 2};
    // Put constraints on u, x
  sfloat Acstr_input_data[2*NINPUTS*NINPUTS] = {0};         // A1*u <= b1
  sfloat Acstr_state_data[2*NSTATES*NSTATES] = {0};         // A2*x <= b2
  // [u_max, -u_min]
  sfloat bcstr_input_data[2*NINPUTS] = {5, 5, 5, 5};  
  // [x_max, -x_min]
  sfloat bcstr_state_data[2*NSTATES] = {6, 8, 3, 2, 
                                        2, 2, 2, 2};
  sfloat input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0};
  sfloat state_dual_data[2 * NSTATES * (NHORIZON)] = {0};
  sfloat goal_dual_data[NSTATES] = {0};

  tiny_LtiModel model;
  tiny_InitLtiModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  Matrix xg = slap_MatrixFromArray(NSTATES, 1, xg_data);

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

  sfloat* Xptr = X_data;
  sfloat* Xref_ptr = Xref_data;
  sfloat* Uptr = U_data;
  sfloat* Uref_ptr = Uref_data;
  sfloat* Kptr = K_data;
  sfloat* dptr = d_data;
  sfloat* Pptr = P_data;
  sfloat* pptr = p_data;
  sfloat* udual_ptr = input_dual_data;
  sfloat* xdual_ptr = state_dual_data;

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
  prob.ncstr_inputs = 1;
  prob.ncstr_states = 1;
  prob.ncstr_goal = 1;

  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 100 * 1e-1);

  prob.Acstr_state = slap_MatrixFromArray(2*NSTATES, NSTATES, Acstr_state_data);
  Matrix upper_half = slap_CreateSubMatrix(prob.Acstr_state, 0, 0, NSTATES, NSTATES);
  Matrix lower_half = slap_CreateSubMatrix(prob.Acstr_state, NSTATES, 0,
                                           NSTATES, NSTATES);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);  
  prob.Acstr_input = slap_MatrixFromArray(2*NINPUTS, NINPUTS, Acstr_input_data);
  upper_half = slap_CreateSubMatrix(prob.Acstr_input, 0, 0, NINPUTS, NINPUTS);
  lower_half = slap_CreateSubMatrix(prob.Acstr_input, NINPUTS, 0,
                                           NINPUTS, NINPUTS);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);                                             
  prob.bcstr_state = slap_MatrixFromArray(2*NSTATES, 1, bcstr_state_data);
  prob.bcstr_input = slap_MatrixFromArray(2*NINPUTS, 1, bcstr_input_data);
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

  solver.max_outer_iters = 10;

  int temp_size = 2*NSTATES * (2*NSTATES + 2*NSTATES + 2)
                + (NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1);
  sfloat temp_data[temp_size];  // temporary data, should not be changed

  tiny_MpcLti(X, U, &prob, &solver, model, 0, temp_data);

  // ========== Test ==========
  for (int k = 0; k < NHORIZON - 1; ++k) {
    tiny_Print(U[k]);
    for (int i = 0; i < NSTATES; ++i) {
      TEST(X[k].data[i] < xmax_data[i] + solver.cstr_tol);
      TEST(X[k].data[i] > xmin_data[i] - solver.cstr_tol);
    }
    for (int i = 0; i < NINPUTS; ++i) {
      TEST(U[k].data[i] > umin_data[i] - solver.cstr_tol);
      TEST(U[k].data[i] < umax_data[i] + solver.cstr_tol);
    }
  }
  tiny_Print(X[NHORIZON - 1]);
  TEST(SumOfSquaredError(X[NHORIZON - 1].data, xg_data, NSTATES) <
       solver.cstr_tol);
}

int main() {
  MpcLtiTest();
  PrintTestResult();
  return TestResult();
}
