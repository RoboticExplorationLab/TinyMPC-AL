// MPC
// Scenerio: Drive bicycle to track references with constraints.

// === BETTER TURN OF GOAL_CONSTRAINT IN PROJECT CMAKELISTS.TXT TO PASS ===
// IF BOX CONSTRAINTS OFF, CAN HANDLE GOAL CONSTRAINT
// IF BOX CONSTRAINTS ON, UNLIKELY TO HANDLE GOAL CONSTRAINT
// NO GRADIENT VANISHING/EXPLOSION WHEN NHORIZON = 91 (101 FAILS)
// GREATER NHORIZON, GREATER ITERATION, GREATER CHANCE OF EXPLOSION
// TODO: Let user choose constraints, compile options with #IFDEF

#include "tiny_mpc_ltv.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "tiny_utils.h"
#include "bicycle_5d.h"
#include "data/lqr_ltv_data.h"

#define H 0.1
#define NSTATES 5
#define NINPUTS 2
#define NHORIZON 5
#define NSIM 101

void MpcTest() {
  double x0_data[NSTATES] = {1, -1, 0, 0, 0};
  // double xg_data[NSTATES] = {0};
  // double ug_data[NINPUTS] = {0};
  double Xhrz_data[NSTATES * NHORIZON] = {0};
  double X_data[NSTATES * NSIM] = {0};
  double Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
  double K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  double d_data[NINPUTS * (NHORIZON - 1)] = {0};
  double P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  double p_data[NSTATES * NHORIZON] = {0};
  double A_data[NSTATES * NSTATES * (NHORIZON-1)] = {0};
  double B_data[NSTATES * NINPUTS * (NHORIZON-1)] = {0};
  double f_data[NSTATES * (NHORIZON-1)] = {0};
  double input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0};
  double state_dual_data[2 * NSTATES * (NHORIZON)] = {0};
  double goal_dual_data[NSTATES] = {0};  
  double Q_data[NSTATES * NSTATES] = {0};
  double R_data[NINPUTS * NINPUTS] = {0};
  double Qf_data[NSTATES * NSTATES] = {0};
  double umin_data[NINPUTS] = {-2.1, -1.1};
  double umax_data[NINPUTS] = {2.1, 1.1};
  double xmin_data[NSTATES] = {-100, -100, -100, -4.0, -0.8};
  double xmax_data[NSTATES] = {100, 100, 100, 4.0, 0.8};

  // double umin_data[NINPUTS] = {-5, -2};
  // double umax_data[NINPUTS] = {5, 2};
  // double xmin_data[NSTATES] = {-100, -100, -100, -100, -100};
  // double xmax_data[NSTATES] = {100, 100, 100, 100, 100};
  
  Matrix X[NSIM];
  Matrix Xref[NSIM];
  Matrix Uref[NSIM - 1];
  Matrix Xhrz[NHORIZON];
  Matrix Uhrz[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1]; 
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix A[NHORIZON-1];
  Matrix B[NHORIZON-1];
  Matrix f[NHORIZON-1];
  Matrix input_duals[NHORIZON - 1];
  Matrix state_duals[NHORIZON];

  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  double* Xhrz_ptr = Xhrz_data;
  double* Xptr = X_data;
  double* Xref_ptr = Xref_data;
  double* Uhrz_ptr = Uhrz_data;
  double* Uref_ptr = Uref_data;
  double* Kptr = K_data;
  double* dptr = d_data;
  double* Pptr = P_data;
  double* pptr = p_data;
  double* Aptr = A_data;
  double* Bptr = B_data;
  double* fptr = f_data;
  double* udual_ptr = input_dual_data;
  double* xdual_ptr = state_dual_data;

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
      Aptr += NSTATES*NSTATES;
      B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, Bptr);
      Bptr += NSTATES*NINPUTS;
      f[i]= slap_MatrixFromArray(NSTATES, 1, fptr);
      fptr += NSTATES;
      Uhrz[i] = slap_MatrixFromArray(NINPUTS, 1, Uhrz_ptr);
      // slap_SetConst(U[i], 0.01);
      Uhrz_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
      input_duals[i] = slap_MatrixFromArray(2 * NINPUTS, 1, udual_ptr);
      udual_ptr += 2 * NINPUTS;
    }
    Xhrz[i] = slap_MatrixFromArray(NSTATES, 1, Xhrz_ptr);
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
  model.get_jacobians = tiny_Bicycle5dGetJacobians;  // from Bicycle
  model.A = A;
  model.B = B;
  model.f = f;
  slap_MatrixCopy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 2 * NINPUTS;
  prob.ncstr_states = 2 * NSTATES;
  prob.ncstr_goal = NSTATES;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
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

  solver.max_primal_iters = 10; // Often takes less than 5

  // Absolute formulation
  // At each time step (end earlier as horizon exceeds the end)
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    printf("\n=> k = %d\n", k);

    // === 1. Setup and solve MPC ===

    slap_MatrixCopy(Xhrz[0], X[k]);
    // Update A, B within horizon 
    tiny_UpdateHorizonJacobians(&model, prob);
    // Update reference
    prob.X_ref = &Xref[k];  
    prob.U_ref = &Uref[k];

    // Solve optimization problem using Augmented Lagrange TVLQR
    tiny_MpcLtv(Xhrz, Uhrz, &prob, &solver, model, 0);

    // Want to print out solution?
    for (int i = 0; i < NHORIZON; ++i) {
      // printf("ex[%d] = %.4f\n", i, slap_MatrixNormedDifference(Xhrz[i], Xref[k+i]));
      // tiny_Print(slap_Transpose(X[k]));
    }

    // === 2. Simulate dynamics using the first control solution ===

    tiny_Bicycle5dNonlinearDynamics(&X[k+1], X[k], Uhrz[0]);

    // Test control constraints here (as we didn't save U)
    TEST(slap_NormInf(Uhrz[0]) < slap_NormInf(prob.u_max) + solver.cstr_tol);
  }

  for (int k = 0; k < NSIM-NHORIZON-1; ++k) {
    printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(X[k], Xref[k]));
    // tiny_Print(slap_Transpose(X[k]));
  }

  // ========== Test ==========
  // Test state constraints
  for (int k = 0; k < NSIM-NHORIZON-1; ++k) {
    for (int i = 0; i < NSTATES; ++i) {
      TEST(X[k].data[i] < xmax_data[i] + solver.cstr_tol);
      TEST(X[k].data[i] > xmin_data[i] - solver.cstr_tol);
    }
  }
  // Test tracking performance
  for (int k = NSIM-NHORIZON-5; k < NSIM-NHORIZON; ++k) {
    TEST(slap_MatrixNormedDifference(X[k], Xref[k]) < 0.2);
  }
  // --------------------------
}

int main() {
  MpcTest();
  PrintTestResult();
  return TestResult();
}