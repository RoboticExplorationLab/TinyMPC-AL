//    __  _             __  _______  ______
//   / /_(_)___  __  __/  |/  / __ \/ ____/
//  / __/ / __ \/ / / / /|_/ / /_/ / /     
// / /_/ / / / / /_/ / /  / / ____/ /___   
// \__/_/_/ /_/\__, /_/  /_/_/    \____/   
//            /____/                                 

// This example run MPC in setup() once because it tracks a fixed smooth trajectory.

#include <Wire.h>
#include <elapsedMillis.h>
#include <stdlib.h>

#include "slap_arduino.h"
#include "tinympc_arduino.h"
#include "bicycle_5d.h"
#include "data/lqr_ltv_data.h"

#define H         0.1       // dt
#define NSTATES   5         // no. of states
#define NINPUTS   2         // no. of controls
#define NHORIZON  10        // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM      101       // simulation steps (fixed with reference data)

#define NOISE(percent) (((2 * ((float)rand() / RAND_MAX)) - 1)/100*percent)

uint64_t run_time = 0;      /* For benchmark */
uint64_t max_run_time = 0;
char bufferTxSer[100];      /* For serial printing */

void setup() {
  Serial.begin(115200);
  while(!Serial) {}
  pinMode(LED_BUILTIN, OUTPUT);
  Serial.println("=============");
  Serial.println("Start problem");

  // ===== Created data =====
  sfloat x0_data[NSTATES] = {1, -1, 0.1, 0, 0};  // initial state (off-track)
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

  // Put constraints on u and x
  sfloat Acstr_input_data[2*NINPUTS*NINPUTS] = {0};         // A1*u <= b1
  sfloat Acstr_state_data[2*NSTATES*NSTATES] = {0};         // A2*x <= b2
  // [u_max, -u_min]
  sfloat bcstr_input_data[2*NINPUTS] = {2.2, 0.9, 2.2, 0.9};  
  // [x_max, -x_min]
  sfloat bcstr_state_data[2*NSTATES] = {100, 100, 100, 4.2, 0.7, 
                                        100, 100, 100, 4.2, 0.7};

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
  sfloat* Xref_ptr = Xref_data;  // Xref defined inside data folder
  sfloat* Uhrz_ptr = Uhrz_data;
  sfloat* Uref_ptr = Uref_data;  // Uref defined inside data folder
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
      tiny_Bicycle5dGetJacobians;  // have analytical functions to compute
                                   // Jacobians, or you can assign manually for
                                   // each time step
  model.get_nonlinear_dynamics =
      tiny_Bicycle5dNonlinearDynamics;  // have dynamics

  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;

  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);

  // Set up constraints
  prob.ncstr_inputs = 1;
  prob.ncstr_states = 1;

  prob.Acstr_state =
      slap_MatrixFromArray(2 * NSTATES, NSTATES, Acstr_state_data);
  Matrix upper_half =
      slap_CreateSubMatrix(prob.Acstr_state, 0, 0, NSTATES, NSTATES);
  Matrix lower_half = slap_CreateSubMatrix(prob.Acstr_state, NSTATES, 0,
                                           NSTATES, NSTATES);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  prob.Acstr_input =
      slap_MatrixFromArray(2 * NINPUTS, NINPUTS, Acstr_input_data);
  upper_half =
      slap_CreateSubMatrix(prob.Acstr_input, 0, 0, NINPUTS, NINPUTS);
  lower_half = slap_CreateSubMatrix(prob.Acstr_input, NINPUTS, 0,
                                    NINPUTS, NINPUTS);
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

  solver.max_outer_iters = 2;   // Often takes less than 5, even still work fine with 2
  solver.cstr_tol = 1e-2;       // AL convergence tolerance (on constraints)

  int temp_size = 2 * NSTATES * (2 * NSTATES + 2 * NSTATES + 2) +
                  (NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1);
  sfloat temp_data[temp_size];
memset(temp_data, 0, sizeof(temp_data));  // temporary data, should not be changed
  srand(1);  // random seed

  // ===== Absolute formulation =====
  // Warm-starting since horizon data is reused
  // At each time step (stop earlier as horizon exceeds the end)
  // Put this into loop() when you have a unlimited reference/real-world system
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    sprintf(bufferTxSer, "k = %d", k);
    Serial.println(bufferTxSer);
    // === 1. Setup and solve MPC ===
    X[k].data[0] += X[k].data[0] * NOISE(1);  // noise within 1% of current X
    X[k].data[1] += X[k].data[1] * NOISE(1);
    X[k].data[2] += X[k].data[2] * NOISE(1);
    X[k].data[3] += X[k].data[3] * NOISE(1);
    X[k].data[4] += X[k].data[4] * NOISE(1);
    slap_Copy(Xhrz[0], X[k]);  // update current measurement

    // Update A, B within horizon (as we have Jacobians function)
    tiny_UpdateHorizonJacobians(&model, prob);

    // Update reference
    prob.X_ref = &Xref[k];
    prob.U_ref = &Uref[k];

    run_time = micros();
    
    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_MpcLtv(Xhrz, Uhrz, &prob, &solver, model, 0, temp_data);

    run_time = (micros() - run_time);   
    sprintf(bufferTxSer, "  MPC run time: %.3f (ms)", ((float)run_time)/1000);
    Serial.println(bufferTxSer);
    max_run_time < run_time? max_run_time = run_time: 0;

    // === 2. Simulate dynamics using the first control solution ===
    tiny_Bicycle5dNonlinearDynamics(&X[k + 1], X[k], Uhrz[0]);

    // tiny_ShiftFill(Uhrz, NHORIZON - 1);  // doesn't matter

    // Print norm of tracking error
    sprintf(bufferTxSer, "  ||ex[%d]|| = %.4f", k, slap_NormedDifference(X[k], Xref[k]));
    Serial.println(bufferTxSer);
  }
  sprintf(bufferTxSer, "\nMPC max run time: %.3f (ms)", ((float)max_run_time)/1000);
  Serial.println(bufferTxSer);
  Serial.println("End of problem");
  Serial.println("==============");
}

void loop() {
  // put your main code here, to run repeatedly:

}
