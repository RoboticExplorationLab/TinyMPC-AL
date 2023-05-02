// Test AL-TVLQR
// Scenerio: Drive bicycle to track references with constraints.

// === BETTER TURN OFF GOAL_CONSTRAINT IN PROJECT CMAKELISTS.TXT TO PASS ===
// IF BOX CONSTRAINTS OFF, CAN HANDLE GOAL CONSTRAINT
// IF BOX CONSTRAINTS ON, UNLIKELY TO HANDLE GOAL CONSTRAINT
// NO GRADIENT VANISHING/EXPLOSION WHEN NHORIZON = 65 (MORE MAY FAIL)
// GREATER NHORIZON, GREATER ITERATION, GREATER CHANCE OF EXPLOSION
// TODO: Let user choose constraints, compile options with #IFDEF

#include "bicycle_5d.h"
#include "data/lqr_ltv_data.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/mpc_ltv.h"
#include "tinympc/utils.h"

#define H 0.1
#define NSTATES 5
#define NINPUTS 2
#define NHORIZON 60

sfloat x0_data[NSTATES] = {1, -1, 0, 0, 0};
sfloat xg_data[NSTATES] = {0};
sfloat ug_data[NINPUTS] = {0};
sfloat Q_data[NSTATES * NSTATES] = {0};
sfloat R_data[NINPUTS * NINPUTS] = {0};
sfloat Qf_data[NSTATES * NSTATES] = {0};
// Put constraints on u, x4, x5
sfloat Acstr_input_data[2 * NINPUTS * NINPUTS] = {0};
sfloat Acstr_state_data[2 * NSTATES * NSTATES] = {0};
// [u_max, -u_min]
sfloat bcstr_input_data[2 * NINPUTS] = {2.0, 0.9, 2.0, 0.9};
// [x_max, -x_min]
sfloat bcstr_state_data[2 * NSTATES] = {100, 100, 100, 4.0, 0.55,
                                        100, 100, 100, 4.0, 0.55};

// sfloat umin_data[NINPUTS] = {-5, -2};
// sfloat umax_data[NINPUTS] = {5, 2};
// sfloat xmin_data[NSTATES] = {-100, -100, -100, -100, -100};
// sfloat xmax_data[NSTATES] = {100, 100, 100, 100, 100};

Matrix X[NHORIZON];
Matrix U[NHORIZON - 1];
Matrix Xref[NHORIZON];
Matrix Uref[NHORIZON - 1];
Matrix K[NHORIZON - 1];
Matrix d[NHORIZON - 1];
Matrix P[NHORIZON];
Matrix p[NHORIZON];
Matrix A[NHORIZON - 1];
Matrix B[NHORIZON - 1];
Matrix f[NHORIZON - 1];
Matrix YU[NHORIZON - 1];
Matrix YX[NHORIZON];

void AbsLqrLtvTest() {
  sfloat X_data[NSTATES * NHORIZON] = {0};
  sfloat U_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  sfloat d_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  sfloat p_data[NSTATES * NHORIZON] = {0};
  sfloat A_data[NSTATES * NSTATES * (NHORIZON - 1)] = {0};
  sfloat B_data[NSTATES * NINPUTS * (NHORIZON - 1)] = {0};
  sfloat f_data[NSTATES * (NHORIZON - 1)] = {0};
  sfloat input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0};
  sfloat state_dual_data[2 * NSTATES * (NHORIZON)] = {0};
  sfloat goal_dual_data[NSTATES] = {0};

  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Settings solver;
  tiny_InitSettings(&solver);

  sfloat* Xptr = X_data;
  sfloat* Xref_ptr = Xref_data;
  sfloat* Uptr = U_data;
  sfloat* Uref_ptr = Uref_data;
  sfloat* Kptr = K_data;
  sfloat* dptr = d_data;
  sfloat* Pptr = P_data;
  sfloat* pptr = p_data;
  sfloat* Aptr = A_data;
  sfloat* Bptr = B_data;
  sfloat* fptr = f_data;
  sfloat* udual_ptr = input_dual_data;
  sfloat* xdual_ptr = state_dual_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      A[i] = slap_MatrixFromArray(NSTATES, NSTATES, Aptr);
      Aptr += NSTATES * NSTATES;
      B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, Bptr);
      Bptr += NSTATES * NINPUTS;
      f[i] = slap_MatrixFromArray(NSTATES, 1, fptr);
      fptr += NSTATES;
      U[i] = slap_MatrixFromArray(NINPUTS, 1, Uptr);
      // slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
      YU[i] = slap_MatrixFromArray(2 * NINPUTS, 1, udual_ptr);
      udual_ptr += 2 * NINPUTS;
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    Xref_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
    YX[i] = slap_MatrixFromArray(2 * NSTATES, 1, xdual_ptr);
    xdual_ptr += 2 * NSTATES;
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  model.get_jacobians = tiny_Bicycle5dGetJacobians;  // from Bicycle
  model.get_nonl_model = tiny_Bicycle5dNonlinearDynamics;
  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 1;
  prob.ncstr_states = 1;
  prob.ncstr_goal = 0;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);
  prob.Acx =
      slap_MatrixFromArray(2 * NSTATES, NSTATES, Acstr_state_data);
  Matrix upper_half =
      slap_CreateSubMatrix(prob.Acx, 0, 0, NSTATES, NSTATES);
  Matrix lower_half =
      slap_CreateSubMatrix(prob.Acx, NSTATES, 0, NSTATES, NSTATES);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  prob.Acu =
      slap_MatrixFromArray(2 * NINPUTS, NINPUTS, Acstr_input_data);
  upper_half = slap_CreateSubMatrix(prob.Acu, 0, 0, NINPUTS, NINPUTS);
  lower_half =
      slap_CreateSubMatrix(prob.Acu, NINPUTS, 0, NINPUTS, NINPUTS);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  prob.bcx = slap_MatrixFromArray(2 * NSTATES, 1, bcstr_state_data);
  prob.bcu = slap_MatrixFromArray(2 * NINPUTS, 1, bcstr_input_data);
  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;
  prob.YU = YU;
  prob.YX = YX;
  prob.YG = slap_MatrixFromArray(NSTATES, 1, goal_dual_data);

  // Absolute formulation
  // Compute and store A, B before solving
  tiny_UpdateHorizonJacobians(&model, prob);

  solver.cstr_tol = 1e-3;
  solver.max_outer_iters = 10;
  int temp_size = 2 * NSTATES * (2 * NSTATES + 2 * NSTATES + 2) +
                  (NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1);
  sfloat temp_data[temp_size];

  memset(temp_data, 0, sizeof(temp_data));
  tiny_MpcLtv(X, U, &prob, &solver, model, 0, temp_data);

  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    // tiny_NonlinearDynamics(&X[k+1], X[k], Uref[k]);
    // PrintMatrix(slap_Transpose(Xref[k]));
    // PrintMatrix(model.B[k]);
  }
  // ========== Test ==========
  for (int k = 0; k < NHORIZON - 1; ++k) {
    // PrintMatrix(U[k]);
    for (int i = 0; i < NSTATES; ++i) {
      TEST(X[k].data[i] < bcstr_state_data[i] + solver.cstr_tol);
      TEST(X[k].data[i] > -bcstr_state_data[i] - solver.cstr_tol);
    }
    for (int i = 0; i < NINPUTS; ++i) {
      TEST(U[k].data[i] > -bcstr_input_data[i] - solver.cstr_tol);
      TEST(U[k].data[i] < bcstr_input_data[i] + solver.cstr_tol);
    }
  }
  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    TEST(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES) < 0.1);
  }
  // --------------------------
}

int main() {
  printf("=== AL LQR LTV Test ===\n");
  AbsLqrLtvTest();
  PrintTestResult();
  return TestResult();
}
