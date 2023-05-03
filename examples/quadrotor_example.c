// MPC
// Scenerio: drive the Crazyflie quadrotor from a state to the origin
//

#include <stdlib.h>

#include "quadrotor.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "time.h"
#include "tinympc/tinympc.h"

// Macro variables
#define H 0.02       // dt
#define NSTATES 12   // no. of states (error state)
#define NINPUTS 4    // no. of controls
#define NHORIZON 15  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 100     // simulation steps (fixed with reference data)

int main() {
  // ===== Created data =====
  sfloat x0_data[NSTATES] = {-0.5, 0.5, -0.5, 0.1, 0, 0,
                             0,    0,   0,    0,   0, 0};  // initial state
  sfloat xg_data[NSTATES] = {0.0};  // goal state if needed
  // sfloat ug_data[NINPUTS] = {0};   // goal input if needed
  sfloat ug_data[NINPUTS] = {0, 0, 0, 0.0};    // goal input if needed
  sfloat Xhrz_data[NSTATES * NHORIZON] = {0};  // save X for one horizon
  sfloat X_data[NSTATES * NSIM] = {0};         // save X for the whole run
  sfloat Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};  // feedback gain
  sfloat d_data[NINPUTS * (NHORIZON - 1)] = {0};            // feedforward gain
  sfloat P_data[NSTATES * NSTATES * (NHORIZON)] = {0};      // cost-to-go func
  sfloat p_data[NSTATES * NHORIZON] = {0};                  // cost-to-go func
  sfloat A_data[NSTATES * NSTATES] = {
      1.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 1.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  1.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, -0.003924f, 0.000000f, 1.000000f, 0.000000f, 0.000000f,
      0.000000f, -0.392400f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.003924f, 0.000000f,  0.000000f, 0.000000f, 1.000000f, 0.000000f,
      0.392400f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 1.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.020000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      1.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.020000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 1.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.020000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  1.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, -0.000013f, 0.000000f, 0.010000f, 0.000000f, 0.000000f,
      0.000000f, -0.001962f, 0.000000f, 1.000000f, 0.000000f, 0.000000f,
      0.000013f, 0.000000f,  0.000000f, 0.000000f, 0.010000f, 0.000000f,
      0.001962f, 0.000000f,  0.000000f, 0.000000f, 1.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.010000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 1.000000f,
  };
  sfloat B_data[NSTATES * NINPUTS] = {
      -0.000019f, -0.000001f, 0.000981f, 0.001264f,  -0.029414f, 0.004771f,
      -0.003847f, -0.000165f, 0.098100f, 0.252748f,  -5.882783f, 0.954290f,
      -0.000001f, -0.000019f, 0.000981f, 0.029044f,  -0.001057f, -0.003644f,
      -0.000138f, -0.003799f, 0.098100f, 5.808852f,  -0.211410f, -0.728857f,
      0.000019f,  0.000001f,  0.000981f, -0.001493f, 0.028771f,  0.001265f,
      0.003763f,  0.000195f,  0.098100f, -0.298680f, 5.754175f,  0.252942f,
      0.000001f,  0.000019f,  0.000981f, -0.028815f, 0.001700f,  -0.002392f,
      0.000222f,  0.003769f,  0.098100f, -5.762921f, 0.340018f,  -0.478376f,
  };
  sfloat f_data[NSTATES] = {0};                                // f in model
  sfloat input_dual_data[2 * NINPUTS * (NHORIZON - 1)] = {0};  // dual vars
  sfloat state_dual_data[2 * NSTATES * (NHORIZON)] = {0};      // dual vars
  sfloat goal_dual_data[NSTATES] = {0};                        // dual vars
  sfloat Q_data[NSTATES * NSTATES] = {0};   // Q matrix in obj
  sfloat R_data[NINPUTS * NINPUTS] = {0};   // R matrix in obj
  sfloat Qf_data[NSTATES * NSTATES] = {0};  // Qf matrix in obj
  sfloat q_data[NSTATES*(NHORIZON-1)] = {0};
  sfloat r_data[NINPUTS*(NHORIZON-1)] = {0};
  sfloat qf_data[NSTATES] = {0};  

  // Put constraints on u
  sfloat Acstr_input_data[2 * NINPUTS * NINPUTS] = {0};  // A1*u <= b1
  sfloat Acstr_state_data[2 * NSTATES * NSTATES] = {0};  // A2*x <= b2
  // [u_max, -u_min]
  sfloat bcstr_input_data[2 * NINPUTS] = {0};
  // [x_max, -x_min]
  sfloat bcstr_state_data[2 * NSTATES] = {0};

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
  Matrix YU[NHORIZON - 1];
  Matrix YX[NHORIZON];
  Matrix q[NHORIZON-1];
  Matrix r[NHORIZON-1];
  Matrix A;
  Matrix B;
  Matrix f;

  // ===== Created tinyMPC struct =====
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
  // tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
  tiny_Settings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all

  tiny_Data data;
  tiny_Info info;
  tiny_Solution soln;
  tiny_Workspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  sfloat temp_data[work.data_size];
  INIT_ZEROS(temp_data);

  tiny_InitTempData(&work, temp_data);

  tiny_InitModelDataArray(&model, &A, &B, &f, A_data, B_data, f_data);

  // ===== Fill in the struct =====
  for (int i = 0; i < NSIM; ++i) {
    if (i < NSIM - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, ug_data);
      // PrintMatrix(Uref[i]);
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, &X_data[i * NSTATES]);
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, xg_data);
    // PrintMatrix(Xref[i]);
  }
  
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      Uhrz[i] = slap_MatrixFromArray(NINPUTS, 1, &Uhrz_data[i * NINPUTS]);
      slap_Copy(Uhrz[i], Uref[i]);  // Initialize U
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES,
                                  &K_data[i * NINPUTS * NSTATES]);
      d[i] = slap_MatrixFromArray(NINPUTS, 1, &d_data[i * NINPUTS]);
      YU[i] = slap_MatrixFromArray(2 * NINPUTS, 1,
                                            &input_dual_data[i * 2 * NINPUTS]);
      q[i] = slap_MatrixFromArray(NSTATES, 1, &q_data[i * NSTATES]);
      r[i] = slap_MatrixFromArray(NINPUTS, 1, &r_data[i * NINPUTS]);
    }
    Xhrz[i] = slap_MatrixFromArray(NSTATES, 1, &Xhrz_data[i * NSTATES]);
    slap_Copy(Xhrz[i], Xref[i]);  // Initialize U
    P[i] =
        slap_MatrixFromArray(NSTATES, NSTATES, &P_data[i * NSTATES * NSTATES]);
    p[i] = slap_MatrixFromArray(NSTATES, 1, &p_data[i * NSTATES]);
    YX[i] =
        slap_MatrixFromArray(2 * NSTATES, 1, &state_dual_data[2 * NSTATES]);
  }

  soln.X = Xhrz;
  soln.U = Uhrz;
  soln.YU = YU;
  soln.YX = YX;
  soln.YG = slap_MatrixFromArray(NSTATES, 1, goal_dual_data);
  soln.K = K;
  soln.d = d;
  soln.P = P;
  soln.p = p;

  data.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);  // check if possible
  data.X_ref = Xref;
  data.U_ref = Uref;
  data.q = q;
  data.r = r;
  data.qf = slap_MatrixFromArray(NSTATES, 1, qf_data);    

  data.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  // slap_SetIdentity(data.Q, 1000e-1);
  sfloat Qdiag[NSTATES] = {10, 10, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1.0};
  slap_SetDiagonal(data.Q, Qdiag, NSTATES);
  data.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(data.R, 0.1);
  data.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_Copy(data.Qf, data.Q);
  slap_ScaleByConst(data.Qf, 1);

  data.Acx =
      slap_MatrixFromArray(2 * NSTATES, NSTATES, Acstr_state_data);
  Matrix upper_half =
      slap_CreateSubMatrix(data.Acx, 0, 0, NSTATES, NSTATES);
  Matrix lower_half =
      slap_CreateSubMatrix(data.Acx, NSTATES, 0, NSTATES, NSTATES);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  data.Acu =
      slap_MatrixFromArray(2 * NINPUTS, NINPUTS, Acstr_input_data);
  upper_half = slap_CreateSubMatrix(data.Acu, 0, 0, NINPUTS, NINPUTS);
  lower_half =
      slap_CreateSubMatrix(data.Acu, NINPUTS, 0, NINPUTS, NINPUTS);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);

  data.bcx = slap_MatrixFromArray(2 * NSTATES, 1, bcstr_state_data);
  slap_SetConst(data.bcx, 100.0);  // x_max = -x_min = 100
  data.bcu = slap_MatrixFromArray(2 * NINPUTS, 1, bcstr_input_data);
  slap_SetConst(data.bcu, 0.5);  // u_max = -u_min = 0.5

  tiny_UpdateLinearCost(&work);

  if (0) {
    printf("\nProblem Info: \n");
    PrintMatrix(work.data->model->A[0]);
    PrintMatrix(work.data->model->B[0]);
    // PrintMatrix(work.data->model->f[0]);
    PrintMatrix(work.data->Q);
    PrintMatrix(work.data->R);
    PrintMatrix(work.data->Qf);
    PrintMatrixT(work.data->x0);
    PrintMatrixT(work.data->X_ref[NHORIZON-5]);
    PrintMatrixT(work.data->U_ref[NHORIZON-5]);
    PrintMatrixT(work.data->q[NHORIZON-5]);
    PrintMatrixT(work.data->r[NHORIZON-5]);
  }

  stgs.en_cstr_goal = 0;
  stgs.en_cstr_inputs = 1;
  stgs.en_cstr_states = 0;
  stgs.max_iter_riccati = 1;
  stgs.max_iter_al = 6;
  stgs.tol_abs_cstr = 1e-2;
  stgs.verbose = 0;
  stgs.reg_min = 1e-6;

  // ===== Absolute formulation =====
  // Warm-starting since horizon data is reused
  // At each time step (stop earlier as horizon exceeds the end)
  slap_Copy(X[0], work.data->x0);  
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    // printf("\n=> k = %d\n", k);
    printf("\n=> ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));

    // === 1. Setup and solve MPC ===
    for (int j = 0; j < NSTATES; ++j) {
      X[k].data[j] += X[k].data[j] * NOISE(0);
    }
    slap_Copy(work.data->x0, X[k]);  // update current measurement

    clock_t start, end;
    double cpu_time_used;
    start = clock();
    // Solve optimization problem using Augmented Lagrangian TVLQR
     tiny_SolveAlLqr(&work);
    // tiny_BackwardPassLti(&prob, solver, model, &Q_temp);
    // tiny_ForwardPassLti(Xhrz, Uhrz, prob, model);
    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    printf("solve time: %f\n", cpu_time_used);

    // Test control constraints here (since we didn't save U)
    // TEST(slap_NormInf(Uhrz[0]) < slap_NormInf(prob.u_max) + solver.cstr_tol);
    PrintMatrixT(Uhrz[0]);
    Matrix pos = slap_CreateSubMatrix(X[k], 0, 0, 3, 1);
    PrintMatrixT(pos);
    // === 2. Simulate dynamics using the first control solution ===
    // tiny_QuadNonlinearDynamics(&X[k + 1], X[k], Uref[k]);
    tiny_QuadNonlinearDynamics(&X[k + 1], X[k], Uhrz[0]);
    // tiny_DynamicsLti(&X[k + 1], X[k], Uref[k], model);
  }

  if (0) {
    for (int k = 0; k < NHORIZON - 1; ++k) {
      // printf("\nk = %d\n", k);
      // PrintMatrixT(Uhrz[k]);
      // PrintMatrixT(Xhrz[k+1]);
      // PrintMatrix(P[k]);
      // PrintMatrixT(Uref[k]);
      // PrintMatrixT(Xref[k+1])
    }
  }

  // ========== Test ==========
  // Test state constraints
  // for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
  //   for (int i = 0; i < NSTATES; ++i) {
  //     TEST(X[k].data[i] < bcstr_state_data[i] + solver.cstr_tol);
  //     TEST(X[k].data[i] > -bcstr_state_data[i] - solver.cstr_tol);
  //   }
  // }
  // // Test tracking performance
  // for (int k = NSIM - NHORIZON - 5; k < NSIM - NHORIZON; ++k) {
  //   TEST(slap_NormedDifference(X[k], Xref[k]) < 0.2);
  // }
  // --------------------------
  PrintTestResult();
  return TestResult();
}
