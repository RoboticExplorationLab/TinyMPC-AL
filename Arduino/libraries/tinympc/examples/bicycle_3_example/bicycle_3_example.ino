//    __  _             __  _______  ______
//   / /_(_)___  __  __/  |/  / __ \/ ____/
//  / __/ / __ \/ / / / /|_/ / /_/ / /     
// / /_/ / / / / /_/ / /  / / ____/ /___   
// \__/_/_/ /_/\__, /_/  /_/_/    \____/   
//            /____/                                 

// This example run MPC in loop() to track a circle reference (for a fixed no. of step).
// Should always use input constraints, otherwise, the yaw angle will be messy (use quaternion?)

#include <Wire.h>                   // in/out pin
#include <elapsedMillis.h>          // timing
#include <stdlib.h>                 // random

#include "slap_arduino.h"           // our linear algebra
#include "tinympc_arduino.h"        // our tinympc
#include "bicycle_3d.h"             // model
#include "data/bicycle3d_track.h"   // data to track

#define H 0.1                       // dt
#define NSTATES 3                   // no. of states
#define NINPUTS 2                   // no. of controls
#define NHORIZON 10                 // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 100                    // simulation steps (fixed with reference data)

#define NOISE(percent) (((2 * ((float)rand() / RAND_MAX)) - 1)/100*percent)   // noise creator with percentage

uint64_t run_time = 0;      /* For benchmark */
uint64_t max_run_time = 0;
char bufferTxSer[100];      /* For serial printing */

// ===== Created data =====
sfloat x0_data[NSTATES] = {-1, -1, 0.2};                    // initial state (off-track)
sfloat Xhrz_data[NSTATES * NHORIZON] = {0};                 // save X for one horizon
sfloat X_data[NSTATES] = {0};                               // save X for one step
sfloat Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
sfloat K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};    // feedback gain
sfloat d_data[NINPUTS * (NHORIZON - 1)] = {0};              // feedforward gain
sfloat P_data[NSTATES * NSTATES * (NHORIZON)] = {0};        // cost-to-go func
sfloat p_data[NSTATES * NHORIZON] = {0};                    // cost-to-go func
sfloat A_data[NSTATES * NSTATES * (NHORIZON - 1)] = {0};    // A in model
sfloat B_data[NSTATES * NINPUTS * (NHORIZON - 1)] = {0};    // B in model
sfloat f_data[NSTATES * (NHORIZON - 1)] = {0};              // f in model
sfloat input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0}; // dual vars
sfloat state_dual_data[2 * NSTATES * (NHORIZON)] = {0};     // dual vars
sfloat Q_data[NSTATES * NSTATES] = {0};                     // Q matrix in obj
sfloat R_data[NINPUTS * NINPUTS] = {0};                     // R matrix in obj
sfloat Qf_data[NSTATES * NSTATES] = {0};                    // Qf matrix in obj

// Put constraints on u and x
sfloat Acstr_input_data[2 * NINPUTS * NINPUTS] = {0};  // A1*u <= b1
sfloat Acstr_state_data[2 * NSTATES * NSTATES] = {0};  // A2*x <= b2
// [u_max, -u_min]
sfloat bcstr_input_data[2 * NINPUTS] = {1.5, 0.6, 1.5, 0.6};
// [x_max, -x_min]
sfloat bcstr_state_data[2 * NSTATES] = {100, 100, 100,
                                        100, 100, 100};

// ===== Created matrices =====
Matrix X;
Matrix Xref[NSIM];        // can be optimized to only save horizon reference
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
tiny_ProblemData prob;
tiny_Solver solver;

int sim_step = 0;
int control_enabled = 1;

void setup() {
  /* serial to display data */
  Serial.begin(115200);
  while(!Serial) {}
  Serial.println("Start setup");
  pinMode(LED_BUILTIN, OUTPUT);

  // ===== Initialize the struct =====  
  tiny_InitLtvModel(&model);
  tiny_InitProblemData(&prob);
  tiny_InitSolver(&solver);

  // ===== Fill in the struct =====
  X = slap_MatrixFromArray(NSTATES, 1, X_data);
  sfloat* Xhrz_ptr = Xhrz_data;
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
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    Xref_ptr += NSTATES;
    // Serial.println(bufferTxSer);
    // sprintf(bufferTxSer, "ref = %.4f\n", Xref[i].data[0]);
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
  model.get_jacobians = tiny_Bicycle3dGetJacobians;  // have analytical functions to compute Jacobians, or you can assign manually for each time step
  model.get_nonlinear_dynamics = tiny_Bicycle3dNonlinearDynamics;  // have full dynamics for simulation
  
  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X, model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;

  // Set up objective function
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);

  // Enable constraints
  prob.ncstr_inputs = 1;
  prob.ncstr_states = 0;

  // Fill in constraints (A in LHS)
  prob.Acstr_state = slap_MatrixFromArray(2*NSTATES, NSTATES, Acstr_state_data);
  Matrix upper_half = slap_CreateSubMatrix(prob.Acstr_state, 0, 0, prob.ninputs, prob.ninputs);
  Matrix lower_half = slap_CreateSubMatrix(prob.Acstr_state, prob.ninputs, 0,
                                           prob.ninputs, prob.ninputs);
  slap_SetIdentity(upper_half, 1);   // Upper half of A is Identity (bound)
  slap_SetIdentity(lower_half, -1);  // Lower half of A is Negative Identity (bound)  
  prob.Acstr_input = slap_MatrixFromArray(2*NINPUTS, NINPUTS, Acstr_input_data);
  upper_half = slap_CreateSubMatrix(prob.Acstr_input, 0, 0, prob.ninputs, prob.ninputs);
  lower_half = slap_CreateSubMatrix(prob.Acstr_input, prob.ninputs, 0,
                                    prob.ninputs, prob.ninputs);
  slap_SetIdentity(upper_half, 1);   // Upper half of A is Identity (bound)
  slap_SetIdentity(lower_half, -1);  // Lower half of A is Negative Identity (bound)  

  // Fill in constraints (b in RHS)
  prob.bcstr_state = slap_MatrixFromArray(2*NSTATES, 1, bcstr_state_data);
  prob.bcstr_input = slap_MatrixFromArray(2*NINPUTS, 1, bcstr_input_data);

  prob.input_duals = input_duals;
  prob.state_duals = state_duals;

  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  solver.max_outer_iters = 10;  // Often takes less than 5
  srand(1);                     // random seed
  Serial.println("Finish setup");
}

void loop() {
  int temp_size = 2*NSTATES * (2*NSTATES + 2*NSTATES + 2)
                  + (NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1);
  sfloat temp_data[temp_size];  // temporary data, should not be changed

  if (control_enabled) { 
    digitalWrite(LED_BUILTIN, HIGH);

    int k = sim_step;   
    sprintf(bufferTxSer, "k = %d", k);
    Serial.println(bufferTxSer);

    // === 1. Setup and solve MPC ===
    X.data[0] += X.data[0] * NOISE(1);  // noise 1% of current X
    X.data[1] += X.data[1] * NOISE(1);
    X.data[2] += X.data[2] * NOISE(1);

    slap_Copy(Xhrz[0], X);              // update current measurement 

    // Update reference
    prob.X_ref = &Xref[k];
    prob.U_ref = &Uref[k];

    // Update A, B within horizon (as we have Jacobians function)
    tiny_UpdateHorizonJacobians(&model, prob);
    
    run_time = micros();
    
    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_MpcLtv(Xhrz, Uhrz, &prob, &solver, model, 0, temp_data);

    run_time = (micros() - run_time);   
    sprintf(bufferTxSer, "  MPC run time: %.3f (ms), max: %.3f (ms)", ((float)run_time)/1000, ((float)max_run_time)/1000);
    Serial.println(bufferTxSer);
    max_run_time < run_time? max_run_time = run_time: 0;

    // === 2. Simulate dynamics using the first control solution ===
    tiny_Bicycle3dNonlinearDynamics(&X, Xhrz[0], Uhrz[0]);

    sprintf(bufferTxSer, "  ||ex|| = %.4f\n", slap_NormedDifference(X, Xref[k+1]));
    Serial.println(bufferTxSer);

    digitalWrite(LED_BUILTIN, LOW);

    sim_step++;

    if (sim_step > NSIM - NHORIZON - 1) {
      control_enabled = 0;
      Serial.println("End of problem");
    }
  }

  digitalWrite(LED_BUILTIN, LOW);
}
