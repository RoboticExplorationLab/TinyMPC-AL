//    __  _             __  _______  ______
//   / /_(_)___  __  __/  |/  / __ \/ ____/
//  / __/ / __ \/ / / / /|_/ / /_/ / /     
// / /_/ / / / / /_/ / /  / / ____/ /___   
// \__/_/_/ /_/\__, /_/  /_/_/    \____/   
//            /____/                                 

// This example run MPC in loop() but will stop after finishing the reference trajectory.
// Better: use Timer/Interupt to call the controller

#include <Wire.h>
#include <elapsedMillis.h>
#include <stdlib.h>

#include "slap_arduino.h"
#include "tinympc_arduino.h"
#include "pointmass.h"

#define H           0.1         // dt
#define NSTATES     2           // no. of states
#define NINPUTS     2           // no. of controls
#define NHORIZON    15          // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM        100         // simulation steps (fixed with reference data)

#define NOISE(percent) (((2 * ((float)rand() / RAND_MAX)) - 1)/100*percent)

uint64_t run_time = 0;      /* For benchmark */
uint64_t max_run_time = 0;
char bufferTxSer[100];      /* For serial printing */
int control_enabled = 1;
int sim_step = 0;

// ===== Created data =====
sfloat x0_data[NSTATES] = {1.0, -2.0};  // initial state
// sfloat x0_data[NSTATES] = {-1, -1, 0.2};  // initial state
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
sfloat Xg_data[NSTATES] = {0.0, 0.0};             // goal state
sfloat Ug_data[NINPUTS] = {0, 0.0};                   // goal control
sfloat Xzeros[NSTATES] = {0};                           // zero state

// Put constraints on u, x
sfloat Acstr_input_data[2 * NINPUTS * NINPUTS] = {0};  // A1*u <= b1
sfloat Acstr_state_data[2 * NSTATES * NSTATES] = {0};  // A2*x <= b2
// [u_max, -u_min]
sfloat bcstr_input_data[2 * NINPUTS] = {0.5, 0.5, 0.5, 0.5};
// [x_max, -x_min]
sfloat bcstr_state_data[2 * NSTATES] = {100, 100,
                                        100, 100.0};
// an obstacle
sfloat x_obs_data[NSTATES] = {0, 0.0};  // obstacle origin (circle)
sfloat r_obs = 1.0;  // obstacle radius 
sfloat r_obs_buff = 0.1;  // since mpc cannot solve till convergence, make some safety buffer
sfloat vecXC_data[NSTATES] = {0.0};
sfloat vecXI_data[NSTATES] = {0.0};

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
Matrix Xg;
Matrix Ug;
Matrix x_obs;

// ===== Created tinyMPC struct =====
tiny_LtvModel model;
tiny_ProblemData prob;
tiny_Solver solver;

void setup() {
  Serial.begin(115200);
  while(!Serial) {}
  pinMode(LED_BUILTIN, OUTPUT);
  Serial.println("=============");
  Serial.println("Start problem");

  tiny_InitLtvModel(&model);
  tiny_InitProblemData(&prob);
  tiny_InitSolver(&solver);

  // ===== Fill in the struct =====
  for (int i = 0; i < NSIM; ++i) {
    if (i < NSIM - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Ug_data);
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, &X_data[i * NSTATES]);
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xg_data);
    // tiny_Print(Xref[i]);
  }
  Xg = slap_MatrixFromArray(NSTATES, 1, Xg_data);
  Ug = slap_MatrixFromArray(NINPUTS, 1, Ug_data);
  x_obs = slap_MatrixFromArray(NSTATES, 1, x_obs_data);
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      A[i] = slap_MatrixFromArray(NSTATES, NSTATES, &A_data[i * NSTATES * NSTATES]);
      B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, &B_data[i * NSTATES * NINPUTS]);
      f[i] = slap_MatrixFromArray(NSTATES, 1, &f_data[i* NSTATES]);
      Uhrz[i] = slap_MatrixFromArray(NINPUTS, 1, &Uhrz_data[i * NINPUTS]);
      slap_Copy(Uhrz[i], Uref[i]);  // Initialize U
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, &K_data[i * NINPUTS * NSTATES]);
      d[i] = slap_MatrixFromArray(NINPUTS, 1, &d_data[i * NINPUTS]);
      input_duals[i] = slap_MatrixFromArray(2 * NINPUTS, 1, &input_dual_data[i * 2 * NINPUTS]);
    }
    Xhrz[i] = slap_MatrixFromArray(NSTATES, 1, &Xhrz_data[i * NSTATES]);
    slap_Copy(Xhrz[i], Xref[i]);  // Initialize U
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, &P_data[i * NSTATES * NSTATES]);
    p[i] = slap_MatrixFromArray(NSTATES, 1, &p_data[i * NSTATES]);
    state_duals[i] = slap_MatrixFromArray(2 * NSTATES, 1, &state_dual_data[2 * NSTATES]);
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);  // delta x0
  model.get_jacobians =
      tiny_PointMassGetJacobians;  // have analytical functions to compute
                                   // Jacobians, or you can assign manually for
                                   // each time step
  model.get_nonlinear_dynamics =
      tiny_PointMassDynamics;  // have dynamics

  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;

  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e0);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e0);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e0);

  // Set up constraints
  prob.ncstr_inputs = 1;
  prob.ncstr_states = 0;
  prob.ncstr_goal = 0;

  prob.Acstr_input = slap_MatrixFromArray(2 * NINPUTS, NINPUTS, Acstr_input_data);
  Matrix upper_half = slap_CreateSubMatrix(prob.Acstr_input, 0, 0, NINPUTS, NINPUTS);
  Matrix lower_half = slap_CreateSubMatrix(prob.Acstr_input, NINPUTS, 0, NINPUTS, NINPUTS);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  prob.bcstr_input = slap_MatrixFromArray(2 * NINPUTS, 1, bcstr_input_data);

  r_obs = r_obs + r_obs_buff;

  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;
  prob.input_duals = input_duals;
  prob.state_duals = state_duals;

  solver.max_outer_iters = 15;   
  solver.cstr_tol = 1e-2;       // AL convergence tolerance (on constraints)

  srand(3);  // random seed
}

void loop() {
  int temp_size = 2 * NSTATES * (2 * NSTATES + 2 * NSTATES + 2) +
                  (NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1);
  sfloat temp_data[temp_size];
  memset(temp_data, 0, sizeof(temp_data));  // temporary data, should not be changed
  
  // ===== Delta formulation =====
  // Warm-starting since horizon data is reused
  // At each time step (stop earlier as horizon exceeds the end)
  // Put this into loop() when you have a unlimited reference/real-world system
  if (control_enabled) {
    int k = sim_step;
    sprintf(bufferTxSer, "k = %d", k);
    Serial.println(bufferTxSer);
    // === 1. Setup and solve MPC ===
    // Update reference
    prob.X_ref = &Xref[k];  // zeros
    prob.U_ref = &Uref[k];  // zeros
    if (k < 100) {
      Xg_data[0] = 0.0;
      Xg_data[1] = 2.0;
    }
    // if (k >= 100) {
    //   Xg_data[0] = -1.0;
    //   Xg_data[1] = -2.0;
    // }
    
    // Update state constraint
    sfloat distance = slap_NormedDifference(X[k], x_obs);
    sprintf(bufferTxSer, "  d2o = %.4f", distance - (r_obs-r_obs_buff));
    Serial.println(bufferTxSer);
    if (distance < 2.0) {
      distance = slap_NormedDifference(Xhrz[5], x_obs);  // number 10 helps
      prob.Acstr_state = slap_MatrixFromArray(2 * NSTATES, NSTATES, Acstr_state_data);
      prob.bcstr_state = slap_MatrixFromArray(2 * NSTATES, 1, bcstr_state_data);
      Matrix upper_half = slap_CreateSubMatrix(prob.Acstr_state, 0, 0, NSTATES, NSTATES);              
      Matrix lower_half = slap_CreateSubMatrix(prob.Acstr_state, NSTATES, 0, NSTATES, NSTATES);                      
      slap_SetIdentity(upper_half, 1);
      slap_SetIdentity(lower_half, -1);
      
      vecXI_data[0] = (x_obs_data[0] - Xhrz[0].data[0]) * (distance - r_obs) / distance;
      vecXI_data[1] = (x_obs_data[1] - Xhrz[0].data[1]) * (distance - r_obs) / distance;
      Acstr_state_data[0] = -2 * (Xhrz[0].data[0] + vecXI_data[0] - x_obs_data[0]);
      Acstr_state_data[1] = -2 * (Xhrz[0].data[1] + vecXI_data[1] - x_obs_data[1]);
      bcstr_state_data[0] = Acstr_state_data[0] * vecXI_data[0] + Acstr_state_data[1] * vecXI_data[1];
      prob.ncstr_states = 1;
    }
    else prob.ncstr_states = 0;

    // Update current measurement, assuming some noises
    // X[k].data[0] += X[k].data[0] * NOISE(2);  // noise within 2% of current X
    // X[k].data[1] += X[k].data[1] * NOISE(2);
    // X[k].data[2] += X[k].data[2] * NOISE(2);
    slap_Copy(Xhrz[0], X[k]);  // update current measurement
    
    // Update A, B within horizon (as we have Jacobians function)
    tiny_UpdateHorizonJacobians(&model, prob);  // since we are using LTI so no needs to update
    run_time = micros();
    
    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_MpcLtv(Xhrz, Uhrz, &prob, &solver, model, 0, temp_data);

    run_time = (micros() - run_time);   
    // sprintf(bufferTxSer, "  MPC run time: %.3f (ms)", ((float)run_time)/1000);
    // Serial.println(bufferTxSer);
    max_run_time < run_time? max_run_time = run_time: 0;

    // === 2. Simulate dynamics using the first control solution ===
    if (Uhrz[0].data[0] >= bcstr_input_data[0]) Uhrz[0].data[0] = bcstr_input_data[0];
    if (Uhrz[0].data[1] >= bcstr_input_data[1]) Uhrz[0].data[1] = bcstr_input_data[1];
    if (Uhrz[0].data[0] <= -bcstr_input_data[0]) Uhrz[0].data[0] = -bcstr_input_data[0];
    if (Uhrz[0].data[1] <= -bcstr_input_data[1]) Uhrz[0].data[1] = -bcstr_input_data[1];
    tiny_PointMassDynamics(&X[k + 1], X[k], Uhrz[0]);
    // tiny_PointMassNonlinearDynamics(&X[k + 1], X[k], Uref[k]);

    // tiny_ShiftFill(Uhrz, NHORIZON - 1);  // doesn't matter

    // Print norm of tracking error
    sprintf(bufferTxSer, "  ||ex[%d]|| = %.4f", k, slap_NormedDifference(X[k], Xg));
    Serial.println(bufferTxSer);
    sprintf(bufferTxSer, "  x = %.4f, %.4f", X[k].data[0], X[k].data[1]);
    Serial.println(bufferTxSer);
    sprintf(bufferTxSer, "  u = %.4f, %.4f", Uhrz[0].data[0], Uhrz[0].data[1]);
    Serial.println(bufferTxSer);
    // tiny_Print(slap_Transpose(Uhrz[0]));  

    sim_step++;
    if (sim_step > NSIM - NHORIZON - 1) {
      control_enabled = 0;
      sprintf(bufferTxSer, "\nMPC max run time: %.3f (ms)", ((float)max_run_time)/1000);
      Serial.println(bufferTxSer);
      Serial.println("End of problem");
      Serial.println("==============");
    }
  }
}
