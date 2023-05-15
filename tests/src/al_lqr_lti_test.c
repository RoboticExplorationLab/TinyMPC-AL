// Task: Test AL-LQR on double integrator with input/state box constraints and
// goal constraint. Scenerio: drive from initial state to goal state.

#include <time.h>

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/al_lqr.h"
#include "tinympc/auxil.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51

void MpcLtiTest() {
  // sfloat tol = 1e-4;
  sfloat A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                      0.1, 0, 1, 0, 0, 0.1, 0, 1};
  sfloat B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
  sfloat f_data[NSTATES] = {0};
  // sfloat x0_data[NSTATES] = {5, 7, 2, -1.4};
  sfloat x0_data[NSTATES] = {1, 0, 0, 1.0};  
  sfloat xg_data[NSTATES] = {0};
  // sfloat xg_data[NSTATES] = {2, 5, -1, 1};
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
  sfloat q_data[NSTATES*(NHORIZON-1)] = {0};
  sfloat r_data[NINPUTS*(NHORIZON-1)] = {0};
  sfloat qf_data[NSTATES] = {0};  
  sfloat umin_data[NINPUTS] = {-5, -5};
  sfloat umax_data[NINPUTS] = {5, 5};
  sfloat xmin_data[NSTATES] = {-2, -2, -2, -2};
  sfloat xmax_data[NSTATES] = {6, 8, 3, 2};
  // Put constraints on u, x
  sfloat Acu_data[2 * NINPUTS * NINPUTS] = {0};  // A1*u <= b1
  sfloat Acx_data[2 * NSTATES * NSTATES] = {0};  // A2*x <= b2
  // [u_max, -u_min]
  sfloat bcu_data[2 * NINPUTS] = {5, 5, 5, 5};
  // [x_max, -x_min]
  sfloat bcx_data[2 * NSTATES] = {6, 8, 3, 2, 2, 2, 2, 2};
  sfloat YU_data[2 * NINPUTS * (NHORIZON - 1)] = {0};
  sfloat YX_data[2 * NSTATES * (NHORIZON)] = {0};
  sfloat YG_data[NSTATES] = {0};

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
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
  Matrix xg = slap_MatrixFromArray(NSTATES, 1, xg_data);

  sfloat* Xref_ptr = Xref_data;
  sfloat* Uref_ptr = Uref_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
    }
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    slap_Copy(Xref[i], xg);
    Xref_ptr += NSTATES;
  }

  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
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

  tiny_InitWorkspaceTempData(&work, temp_data);

  tiny_InitModelFromArray(&model, &A, &B, &f, A_data, B_data, f_data);

  tiny_InitSolnTrajFromArray(&work, X, U, X_data, U_data);
  tiny_InitSolnDualsFromArray(&work, YX, YU, YX_data, YU_data, YG_data);
  tiny_InitSolnGainsFromArray(&work, K, d, P, p, K_data, d_data, P_data, p_data);

  data.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);  
  data.X_ref = Xref;
  data.U_ref = Uref;
  tiny_InitDataQuadCostFromArray(&work, Q_data, R_data, Qf_data);
  slap_SetIdentity(data.Q, 1);  
  slap_SetIdentity(data.R, 1);  
  slap_SetIdentity(data.Qf, 10);
  tiny_InitDataLinearCostFromArray(&work, q, r, q_data, r_data, qf_data);

  tiny_SetInputBound(&work, Acu_data, bcu_data);
  tiny_SetStateBound(&work, Acx_data, bcx_data);

  tiny_UpdateLinearCost(&work);

  if (0) {
    printf("\nProblem Info: \n");
    PrintMatrix(work.data->model->A[0]);
    PrintMatrix(work.data->model->B[0]);
    PrintMatrix(work.data->model->f[0]);
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

  clock_t start, end;
  double cpu_time_used;
  start = clock();
  tiny_SolveAlLqr(&work);
  end = clock();
  cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;
  printf("time: %f\n", cpu_time_used);

  if (0) {
    for (int k = 0; k < NHORIZON - 1; ++k) {
      // printf("\n=>k = %d\n", k);
      // PrintMatrix(p[k]);
      // PrintMatrixT(Xref[k]);
      // PrintMatrixT(U[k]);
      PrintMatrixT(X[k]);
    }
    PrintMatrixT(X[NHORIZON - 1]);
  }  

  // ========== Test ==========
  for (int k = 0; k < NHORIZON - 1; ++k) {
    for (int i = 0; i < NSTATES; ++i) {
      TEST(X[k].data[i] < xmax_data[i] + stgs.tol_abs_cstr);
      TEST(X[k].data[i] > xmin_data[i] - stgs.tol_abs_cstr);
    }
    for (int i = 0; i < NINPUTS; ++i) {
      TEST(U[k].data[i] > umin_data[i] - stgs.tol_abs_cstr);
      TEST(U[k].data[i] < umax_data[i] + stgs.tol_abs_cstr);
    }
  }
  TEST(SumOfSquaredError(X[NHORIZON - 1].data, xg_data, NSTATES) <
       stgs.tol_abs_cstr);
}

int main() {
  printf("=== AL LQR LTI Test ===\n");
  MpcLtiTest();
  PrintTestResult();
  return TestResult();
}
