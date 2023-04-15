// MPC
// Scenerio: Drive bicycle to track references with constraints.
// Check trajopt at test/al_lqr_test/al_lqr_ltv_test.c

// === BETTER TURN OFF GOAL_CONSTRAINT IN PROJECT CMAKELISTS.TXT TO PASS ===
// IF BOX CONSTRAINTS OFF, CAN HANDLE GOAL CONSTRAINT
// IF BOX CONSTRAINTS ON, UNLIKELY TO HANDLE GOAL CONSTRAINT
// DON"T WORRY ABOUT GRADIENT VANISHING/EXPLOSION SINCE SMALL MPC HORIZON
// GREATER NHORIZON, GREATER ITERATION, GREATER CHANCE OF EXPLOSION
// TODO: Let user choose constraints, compile options with #IFDEF

#include "bicycle_3d.h"
#include "data/bicycle3d_track.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "tinympc/tinympc.h"

#define H 0.1        // dt
#define NSTATES 3    // no. of states
#define NINPUTS 2    // no. of controls
#define NHORIZON 10  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 100   // simulation steps (fixed with reference data)

int main() {
  // ===== Created data =====
  sfloat x0_data[NSTATES] = {-1, -1, 0.2};  // initial state
  sfloat Xhrz_data[NSTATES * NHORIZON] = {0};  // save X for one horizon
  sfloat X_data[NSTATES * NSIM] = {0};         // save X for the whole run
  sfloat Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};  // feedback gain
  sfloat d_data[NINPUTS * (NHORIZON - 1)] = {0};            // feedforward gain
  sfloat P_data[NSTATES * NSTATES * (NHORIZON)] = {0};      // cost-to-go func
  sfloat p_data[NSTATES * NHORIZON] = {0};                  // cost-to-go func
  sfloat A_data[NSTATES * NSTATES * (NHORIZON - 1)] = {0};  // A in model
  sfloat B_data[NSTATES * NINPUTS * (NHORIZON - 1)] = {0};  // B in model
  sfloat f_data[NSTATES * (NHORIZON - 1)] = {0};            // f in model
  sfloat input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0};  // dual vars
  sfloat state_dual_data[2 * NSTATES * (NHORIZON)] = {0};      // dual vars
  sfloat Q_data[NSTATES * NSTATES] = {0};   // Q matrix in obj
  sfloat R_data[NINPUTS * NINPUTS] = {0};   // R matrix in obj
  sfloat Qf_data[NSTATES * NSTATES] = {0};  // Qf matrix in obj

  // Put constraints on u, x4, x5
  sfloat Acstr_input_data[2 * NINPUTS * NINPUTS] = {0};  // A1*u <= b1
  sfloat Acstr_state_data[2 * NSTATES * NSTATES] = {0};  // A2*x <= b2
  // [u_max, -u_min]
  sfloat bcstr_input_data[2 * NINPUTS] = {1.5, 0.6, 1.5, 0.6};
  // [x_max, -x_min]
  sfloat bcstr_state_data[2 * NSTATES] = {100, 100, 100,
                                          100, 100, 100};

  // ===== Created matrices =====
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

  // ===== Created tinyMPC struct =====
  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  // ===== Fill in the struct =====
  sfloat* Xhrz_ptr = Xhrz_data;
  sfloat* Xptr = X_data;
  sfloat* Xref_ptr = X_ref_data;  // Xref defined inside data folder
  sfloat* Uhrz_ptr = Uhrz_data;
  sfloat* Uref_ptr = U_ref_data;  // Uref defined inside data folder
  sfloat* Kptr = K_data;
  sfloat* dptr = d_data;
  sfloat* Pptr = P_data;
  sfloat* pptr = p_data;
  sfloat* Aptr = A_data;
  sfloat* Bptr = B_data;
  sfloat* fptr = f_data;
  sfloat* udual_ptr = input_dual_data;
  sfloat* xdual_ptr = state_dual_data;

  for (int i = 0; i < NSIM; ++i) {
    if (i < NSIM - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    Xref_ptr += NSTATES;
    // tiny_Print(Xref[i]);
  }
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      A[i] = slap_MatrixFromArray(NSTATES, NSTATES, Aptr);
      Aptr += NSTATES * NSTATES;
      B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, Bptr);
      Bptr += NSTATES * NINPUTS;
      f[i] = slap_MatrixFromArray(NSTATES, 1, fptr);
      fptr += NSTATES;
      Uhrz[i] = slap_MatrixFromArray(NINPUTS, 1, Uhrz_ptr);
      slap_Copy(Uhrz[i], Uref[i]);  // Initialize U
      Uhrz_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
      input_duals[i] = slap_MatrixFromArray(2 * NINPUTS, 1, udual_ptr);
      udual_ptr += 2 * NINPUTS;
    }
    Xhrz[i] = slap_MatrixFromArray(NSTATES, 1, Xhrz_ptr);
    slap_Copy(Xhrz[i], Xref[i]);  // Initialize U
    Xhrz_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
    state_duals[i] = slap_MatrixFromArray(2 * NSTATES, 1, xdual_ptr);
    xdual_ptr += 2 * NSTATES;
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  model.get_jacobians =
      tiny_Bicycle3dGetJacobians;  // have analytical functions to compute
                                   // Jacobians, or you can assign manually for
                                   // each time step
  model.get_nonlinear_dynamics =
      tiny_Bicycle3dNonlinearDynamics;  // have dynamics

  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 1;
  prob.ncstr_states = 0;

  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);

  prob.Acstr_state =
      slap_MatrixFromArray(2 * NSTATES, NSTATES, Acstr_state_data);
  Matrix upper_half =
      slap_CreateSubMatrix(prob.Acstr_state, 0, 0, prob.ninputs, prob.ninputs);
  Matrix lower_half = slap_CreateSubMatrix(prob.Acstr_state, prob.ninputs, 0,
                                           prob.ninputs, prob.ninputs);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  prob.Acstr_input =
      slap_MatrixFromArray(2 * NINPUTS, NINPUTS, Acstr_input_data);
  upper_half =
      slap_CreateSubMatrix(prob.Acstr_input, 0, 0, prob.ninputs, prob.ninputs);
  lower_half = slap_CreateSubMatrix(prob.Acstr_input, prob.ninputs, 0,
                                    prob.ninputs, prob.ninputs);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  prob.bcstr_state = slap_MatrixFromArray(2 * NSTATES, 1, bcstr_state_data);
  prob.bcstr_input = slap_MatrixFromArray(2 * NINPUTS, 1, bcstr_input_data);

  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;
  prob.input_duals = input_duals;
  prob.state_duals = state_duals;

  solver.max_outer_iters = 10;  // Often takes less than 5
  int temp_size = 2 * NSTATES * (2 * NSTATES + 2 * NSTATES + 2) +
                  (NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1);
  sfloat temp_data[temp_size];  // temporary data, should not be changed

  // ===== Absolute formulation =====
  // Warm-starting since horizon data is reused
  // At each time step (stop earlier as horizon exceeds the end)
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    printf("\n=> k = %d\n", k);
    // === 1. Setup and solve MPC ===

    slap_Copy(Xhrz[0], X[k]);  // update current measurement

    // Update A, B within horizon (as we have Jacobians function)
    tiny_UpdateHorizonJacobians(&model, prob);

    // Update reference
    prob.X_ref = &Xref[k];
    prob.U_ref = &Uref[k];

    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_MpcLtv(Xhrz, Uhrz, &prob, &solver, model, 0, temp_data);

    // Test control constraints here (since we didn't save U)
    // TEST(slap_NormInf(Uhrz[0]) < slap_NormInf(prob.u_max) + solver.cstr_tol);

    // === 2. Simulate dynamics using the first control solution ===
    tiny_Bicycle3dNonlinearDynamics(&X[k + 1], X[k], Uhrz[0]);
    // tiny_DynamicsLtv(&X[k + 1], X[k], Uhrz[0], model, 0);

    printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
  }

  // ========== Test ==========
  // Test state constraints
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    for (int i = 0; i < NSTATES; ++i) {
      // TEST(X[k].data[i] < bcstr_state_data[i] + solver.cstr_tol);
      // TEST(X[k].data[i] > -bcstr_state_data[i] - solver.cstr_tol);
    }
  }
  // Test tracking performance
  for (int k = NSIM - NHORIZON - 5; k < NSIM - NHORIZON; ++k) {
    // TEST(slap_NormedDifference(X[k], Xref[k]) < 0.1);
  }
  // --------------------------

  PrintTestResult();
  return TestResult();
}
