// MPC
// Scenerio: Drive bicycle to track references with constraints.
// Check trajopt at test/al_lqr_test/al_lqr_ltv_test.c

// === BETTER TURN OFF GOAL_CONSTRAINT IN PROJECT CMAKELISTS.TXT TO PASS ===
// IF BOX CONSTRAINTS OFF, CAN HANDLE GOAL CONSTRAINT
// IF BOX CONSTRAINTS ON, UNLIKELY TO HANDLE GOAL CONSTRAINT
// DON"T WORRY ABOUT GRADIENT VANISHING/EXPLOSION SINCE SMALL MPC HORIZON
// GREATER NHORIZON, GREATER ITERATION, GREATER CHANCE OF EXPLOSION
// TODO: Let user choose constraints, compile options with #IFDEF

#include "planar_quadrotor.h"
#include "data/lqr_ltv_data.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "tinympc/tinympc.h"


#define H 0.05
#define NSTATES 6
#define NINPUTS 2
#define NHORIZON 10
#define NSIM 201

int main() {

  // sfloat g = 9.81; // m/s^2
  // sfloat m = 1.0; // kg

  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  model.nstates = NSTATES;
  model.ninputs = NINPUTS;
  model.dt = H;

  sfloat x0_data[NSTATES] = {1, 3, .5, 0, 0, 0};
  // sfloat xg_data[NSTATES] = {0};
  // sfloat ug_data[NINPUTS] = {0};
  sfloat Xhrz_data[NSTATES * NHORIZON] = {0};
  sfloat X_data[NSTATES * NSIM] = {0};
  sfloat Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
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
  
  sfloat Q_data[NSTATES * NSTATES] = {0};
  sfloat R_data[NINPUTS * NINPUTS] = {0};
  sfloat Qf_data[NSTATES * NSTATES] = {0};

  // Put constraints on u, x4, x5
  // sfloat umin_data[NINPUTS] = {0.2*m*g, 0.2*m*g};
  // sfloat umax_data[NINPUTS] = {0.6*m*g, 0.6*m*g};
  sfloat umin_data[NINPUTS] = {-99999, -99999};
  sfloat umax_data[NINPUTS] = {99999, 99999};
  sfloat xmin_data[NSTATES] = {-99999, 0, -99999, -99999, -99999, -99999};
  sfloat xmax_data[NSTATES] = {99999, 99999, 99999, 99999, 99999, 99999};

  Matrix X[NSIM];

  Matrix Xref[NSIM];
  Matrix Uref[NSIM - 1];
  Matrix Xhrz[NHORIZON];
  Matrix Uhrz[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix A[NHORIZON - 1];
  Matrix B[NHORIZON - 1];
  Matrix f[NHORIZON - 1];
  Matrix input_duals[NHORIZON - 1];
  Matrix state_duals[NHORIZON];

  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  model.get_jacobians = tiny_PQuadGetJacobians;
  model.get_nonlinear_dynamics = tiny_PQuadNonlinearDynamics;
  model.A = A;
  model.B = B;
  model.f = f;

  X[0] = slap_MatrixFromArray(NSTATES, 1, x0_data);

  for (int i = 0; i < NSIM; ++i) {
    if (i < NSIM - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, &Uref_data[i*NINPUTS]);
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, &X_data[i*NSTATES]);
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, &Xref_data[i*NSTATES]);
  }
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      A[i] = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
      B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
      f[i] = slap_MatrixFromArray(NSTATES, 1, f_data);
      Uhrz[i] = slap_MatrixFromArray(NINPUTS, 1, Uhrz_data);
      slap_Copy(Uhrz[i], Uref[i]);  // Initialize U
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, K_data);
      d[i] = slap_MatrixFromArray(NINPUTS, 1, d_data);
      input_duals[i] = slap_MatrixFromArray(2 * NINPUTS, 1, input_dual_data);
    }
    Xhrz[i] = slap_MatrixFromArray(NSTATES, 1, Xhrz_data);
    slap_Copy(Xhrz[i], Xref[i]);  // Initialize U
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, P_data);
    p[i] = slap_MatrixFromArray(NSTATES, 1, p_data);
    state_duals[i] = slap_MatrixFromArray(2 * NSTATES, 1, state_dual_data);
  }

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 2 * NINPUTS;
  prob.ncstr_states = 2 * NSTATES;
  prob.ncstr_goal = 0;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);
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

  solver.max_primal_iters = 10;  // Often takes less than 5

  // Absolute formulation
  // Warm-starting since horizon data is reused
  // At each time step (stop earlier as horizon exceeds the end)
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    // printf("\n=> k = %d\n", k);
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    // === 1. Setup and solve MPC ===

    slap_Copy(Xhrz[0], X[k]);
    // Update A, B within horizon
    tiny_UpdateHorizonJacobians(&model, prob);
    // Update reference
    prob.X_ref = &Xref[k];
    prob.U_ref = &Uref[k];

    // Solve optimization problem using Augmented Lagrangian TVLQR, benchmark
    // this
    tiny_MpcLtv(Xhrz, Uhrz, &prob, &solver, model, 0);

    // Test control constraints here (since we didn't save U)
    TEST(slap_NormInf(Uhrz[0]) < slap_NormInf(prob.u_max) + solver.cstr_tol);

    // === 2. Simulate dynamics using the first control solution ===

    // Clamping control would not effect since our solution is feasible
    tiny_ClampMatrix(&Uhrz[0], prob.u_min, prob.u_max);
    tiny_PQuadNonlinearDynamics(&X[k + 1], X[k], Uhrz[0]);
  }

  // ========== Test ==========
  // Test state constraints
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    for (int i = 0; i < NSTATES; ++i) {
      TEST(X[k].data[i] < xmax_data[i] + solver.cstr_tol);
      TEST(X[k].data[i] > xmin_data[i] - solver.cstr_tol);
    }
  }
  // Test tracking performance
  for (int k = NSIM - NHORIZON - 5; k < NSIM - NHORIZON; ++k) {
    TEST(slap_NormedDifference(X[k], Xref[k]) < 0.1);
  }
  // --------------------------

  PrintTestResult();
  return TestResult();
}
