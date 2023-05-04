#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/auxil.h"
#include "tinympc/model.h"

#define NSTATES 2
#define NINPUTS 1
#define NHORIZON 3

sfloat dt = 0.1;
sfloat A_data[NSTATES * NSTATES] = {1, 0, 1, 1};  // NOLINT
sfloat B_data[NSTATES * NINPUTS] = {1, 2};        // NOLINT
sfloat f_data[NSTATES] = {4, 5};                  // NOLINT
// sfloat t = 1.0;
sfloat Q_data[NSTATES * NSTATES] = {1, 0, 0, 1};  // NOLINT
sfloat R_data[NINPUTS * NINPUTS] = {1};           // NOLINT
sfloat q_data[NSTATES] = {0.1, 0.2};              // NOLINT
sfloat r_data[NINPUTS] = {-0.6};                  // NOLINT
sfloat Qf_data[NSTATES * NSTATES] = {0};
// sfloat u_max_data[NINPUTS] = {1.1};
// sfloat u_min_data[NINPUTS] = {-1.1};
// sfloat x_max_data[NSTATES] = {1.6, 1.7};
// sfloat x_min_data[NSTATES] = {-1.6, -1.7};
sfloat x_ref_data[NSTATES * NHORIZON] = {0.2, 1.1, 2.5, 3.7, 2.1, 4.5};
sfloat u_ref_data[NINPUTS * (NHORIZON - 1)] = {1, 2};
sfloat x0_data[NSTATES] = {0.1, 0.2};
sfloat u0_data[NSTATES] = {1.6};
sfloat Kd_data[NINPUTS * (NSTATES + 1) * (NHORIZON - 1)] = {0};
sfloat Pp_data[NSTATES * (NSTATES + 1) * NHORIZON] = {0};
sfloat reg = 1e-6;
sfloat input_duals_data[2 * NINPUTS * (NHORIZON - 1)] = {1, 2, 3, 4};
sfloat state_duals_data[2 * NSTATES * (NHORIZON)] = {1, 2, 3, 4, 5, 6,
                                                     7, 8, 2, 4, 5, 6};
sfloat goal_duals_data[NSTATES] = {1, 2};
sfloat reg_min = 1;
sfloat reg_max = 100;
sfloat penalty_max = 1e5;
sfloat penalty_mul = 1;
int max_outer_iters = 100;
int max_search_iters = 10;

// void KnotPointTest() {
//   const sfloat tol = 1e-6;
//   tiny_KnotPoint z;
//   tiny_InitKnotPoint(&z);
//   tiny_KnotPoint Z[NHORIZON];
//   z.dt = dt;
//   z.t = t;
//   z.x = slap_MatrixFromArray(NSTATES, 1, x0_data);
//   z.u = slap_MatrixFromArray(NINPUTS, 1, u0_data);

//   TEST(z.dt == dt);
//   TEST(z.t == t);
//   TEST(z.x.rows == NSTATES);
//   TEST(z.x.cols == 1);
//   TEST(SumOfSquaredError(x0_data, z.x.data, z.x.rows) < tol);
//   TEST(z.u.rows == NINPUTS);
//   TEST(z.u.cols == 1);
//   TEST(SumOfSquaredError(u0_data, z.u.data, z.u.rows) < tol);

//   for (int i = 0; i < NHORIZON; ++i) {
//     tiny_InitKnotPoint(&Z[i]);
//     Z[i].dt = dt;
//     Z[i].t = t + i;
//     Z[i].x = slap_MatrixFromArray(NSTATES, 1, x0_data);
//     Z[i].u = slap_MatrixFromArray(NINPUTS, 1, u0_data);
//   }

//   for (int i = 0; i < NHORIZON; ++i) {
//     TEST(Z[i].dt == dt);
//     TEST(Z[i].t == (t + i));
//     TEST(Z[i].x.rows == NSTATES);
//     TEST(Z[i].x.cols == 1);
//     TEST(SumOfSquaredError(Z[i].x.data, z.x.data, z.x.rows) < tol);
//     TEST(Z[i].u.rows == NINPUTS);
//     TEST(Z[i].u.cols == 1);
//     TEST(SumOfSquaredError(Z[i].u.data, z.u.data, z.u.rows) < tol);
//   }
// }

void SettingsTest() {
  tiny_Settings settings;
  tiny_InitSettings(&settings);

  TEST(settings.reg_max == (sfloat)REG_MAX);
  TEST(settings.max_iter_al == MAX_ITER_AL);

  settings.reg_max = (sfloat)10.5;
  TEST(settings.reg_max == (sfloat)10.5);
}

void SolutionTest() {
  int n = NSTATES;
  int m = NINPUTS;
  int N = NHORIZON;
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, dt);
  tiny_Settings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  tiny_Data data;
  tiny_Info info;
  tiny_Solution soln;
  tiny_Workspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  int true_size = n*N + m*(N-1) + (m+1)*n*(N-1) + (n+1)*n*N;
  true_size += m*(N-1)*2 + n*N*2 + n;
  // TEST(work.soln->data_size == true_size);
}

void DataTest() {
  int n = NSTATES;
  int m = NINPUTS;
  int N = NHORIZON;
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, dt);
  tiny_Settings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  tiny_Data data;
  tiny_Info info;
  tiny_Solution soln;
  tiny_Workspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  int true_size = n + n*n*2 + m*m + (N-1)*(n + m) + m + N*n + (N-1)*m;
  true_size += 2*m*m + 2*m + 2*n*n + 2*n;
  // TEST(work.data->data_size == true_size);
}

void WorkspaceTest() {
  int n = NSTATES;
  int m = NINPUTS;
  // int N = NHORIZON;
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, dt);
  tiny_Settings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  tiny_Data data;
  tiny_Info info;
  tiny_Solution soln;
  tiny_Workspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);

  TEST(work.first_run == 1);

  int true_size =  2 * n * (2 * n + 2 * n + 2) + (n + m) * (n + m + 1);
  TEST(work.data_size == true_size);
  
  sfloat temp_data[work.data_size];
  INIT_ZEROS(temp_data);

  tiny_InitWorkspaceTempData(&work, temp_data);

  TEST((work.YU_hat.rows == 2 * m) && (work.YU_hat.cols == 1));
  TEST((work.cu2.rows == 2 * m) && (work.cu2.cols == 1));
}
  
// void ProblemDataTest() {
//   const sfloat tol = 1e-6;
//   Matrix U_ref[NHORIZON - 1];
//   Matrix X_ref[NHORIZON];
//   Matrix K[NHORIZON - 1];
//   Matrix d[NHORIZON - 1];
//   Matrix P[NHORIZON];
//   Matrix p[NHORIZON];
//   Matrix YU[NHORIZON - 1];
//   Matrix YX[NHORIZON];
//   sfloat* uptr = u_ref_data;
//   sfloat* xptr = x_ref_data;
//   sfloat* Kd_ptr = Kd_data;
//   sfloat* Pp_ptr = Pp_data;
//   sfloat* input_dual_ptr = input_duals_data;
//   sfloat* state_dual_ptr = state_duals_data;
//   for (int i = 0; i < NHORIZON; ++i) {
//     if (i < NHORIZON - 1) {
//       U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
//       uptr += NINPUTS;
//       K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd_ptr);
//       Kd_ptr += NINPUTS * NSTATES;
//       d[i] = slap_MatrixFromArray(NINPUTS, 1, Kd_ptr);
//       Kd_ptr += NINPUTS;
//       YU[i] = slap_MatrixFromArray(NINPUTS, 1, input_dual_ptr);
//       input_dual_ptr += NINPUTS;
//     }
//     X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
//     xptr += NSTATES;
//     P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pp_ptr);
//     Pp_ptr += NSTATES * NSTATES;
//     p[i] = slap_MatrixFromArray(NSTATES, 1, Pp_ptr);
//     Pp_ptr += NSTATES;
//     YX[i] = slap_MatrixFromArray(NSTATES, 1, state_dual_ptr);
//     state_dual_ptr += NSTATES;
//   }

//   tiny_ProblemData prob;
//   tiny_InitProblemData(&prob);

//   prob.nstates = NSTATES;
//   prob.ninputs = NINPUTS;
//   prob.nhorizon = NHORIZON;
//   prob.ncstr_goal = NSTATES;
//   prob.ncstr_inputs = 2 * NINPUTS * (NHORIZON - 1);
//   prob.ncstr_states = 2 * NSTATES * NHORIZON;
//   prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
//   prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
//   prob.q = slap_MatrixFromArray(NSTATES, 1, q_data);
//   prob.r = slap_MatrixFromArray(NINPUTS, 1, r_data);
//   prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
//   // prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
//   // prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
//   // prob.x_max = slap_MatrixFromArray(NSTATES, 1, x_max_data);
//   // prob.x_min = slap_MatrixFromArray(NSTATES, 1, x_min_data);
//   prob.X_ref = X_ref;
//   prob.U_ref = U_ref;
//   prob.dt = dt;
//   prob.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
//   prob.K = K;
//   prob.d = d;
//   prob.P = P;
//   prob.p = p;
//   prob.YU = YU;
//   prob.YX = YX;
//   prob.YG = slap_MatrixFromArray(NSTATES, 1, goal_duals_data);

//   TEST(prob.nstates == NSTATES);
//   TEST(prob.ninputs == NINPUTS);
//   TEST(prob.nhorizon == NHORIZON);
//   TEST(prob.ncstr_inputs == 2 * NINPUTS * (NHORIZON - 1));
//   TEST(prob.ncstr_states == 2 * NSTATES * NHORIZON);
//   TEST(prob.dt == dt);
//   TEST(SumOfSquaredError(prob.Q.data, Q_data, NSTATES * NSTATES) < tol);
//   TEST(SumOfSquaredError(prob.R.data, R_data, NINPUTS * NINPUTS) < tol);
//   TEST(SumOfSquaredError(prob.q.data, q_data, NSTATES) < tol);
//   TEST(SumOfSquaredError(prob.r.data, r_data, NINPUTS) < tol);
//   TEST(SumOfSquaredError(prob.Qf.data, Qf_data, NSTATES * NSTATES) < tol);
//   // TEST(SumOfSquaredError(prob.u_max.data, u_max_data, NINPUTS) < tol);
//   // TEST(SumOfSquaredError(prob.u_min.data, u_min_data, NINPUTS) < tol);
//   // TEST(SumOfSquaredError(prob.x_max.data, x_max_data, NSTATES) < tol);
//   // TEST(SumOfSquaredError(prob.x_min.data, x_min_data, NSTATES) < tol);
//   TEST(SumOfSquaredError(prob.x0.data, x0_data, NSTATES) < tol);
//   TEST(SumOfSquaredError(prob.YG.data, goal_duals_data, NSTATES) < tol);
//   uptr = u_ref_data;
//   xptr = x_ref_data;
//   Kd_ptr = Kd_data;
//   Pp_ptr = Pp_data;
//   input_dual_ptr = input_duals_data;
//   state_dual_ptr = state_duals_data;
//   for (int i = 0; i < NHORIZON; ++i) {
//     if (i < NHORIZON - 1) {
//       TEST(prob.U_ref[i].rows == NINPUTS);
//       TEST(prob.U_ref[i].cols == 1);
//       TEST(SumOfSquaredError(prob.U_ref[i].data, uptr, NINPUTS) < tol);
//       uptr += NINPUTS;
//       TEST(prob.K[i].rows == NINPUTS);
//       TEST(prob.K[i].cols == NSTATES);
//       TEST(SumOfSquaredError(prob.K[i].data, Kd_ptr, NINPUTS * NSTATES) < tol);
//       Kd_ptr += NINPUTS * NSTATES;
//       TEST(prob.d[i].rows == NINPUTS);
//       TEST(prob.d[i].cols == 1);
//       TEST(SumOfSquaredError(prob.d[i].data, Kd_ptr, NINPUTS) < tol);
//       Kd_ptr += NINPUTS;
//       TEST(prob.YU[i].rows == NINPUTS);
//       TEST(prob.YU[i].cols == 1);
//       TEST(SumOfSquaredError(prob.YU[i].data, input_dual_ptr,
//                              NINPUTS) < tol);
//       input_dual_ptr += NINPUTS;
//     }
//     TEST(prob.X_ref[i].rows == NSTATES);
//     TEST(prob.X_ref[i].cols == 1);
//     TEST(SumOfSquaredError(prob.X_ref[i].data, xptr, NSTATES) < tol);
//     xptr += NSTATES;
//     TEST(prob.P[i].rows == NSTATES);
//     TEST(prob.P[i].cols == NSTATES);
//     TEST(SumOfSquaredError(prob.P[i].data, Pp_ptr, NSTATES * NSTATES) < tol);
//     Pp_ptr += NSTATES * NSTATES;
//     TEST(prob.p[i].rows == NSTATES);
//     TEST(prob.p[i].cols == 1);
//     TEST(SumOfSquaredError(prob.p[i].data, Pp_ptr, NSTATES) < tol);
//     Pp_ptr += NSTATES;
//     TEST(prob.YX[i].rows == NSTATES);
//     TEST(prob.YX[i].cols == 1);
//     TEST(SumOfSquaredError(prob.YX[i].data, state_dual_ptr, NSTATES) <
//          tol);
//     state_dual_ptr += NSTATES;
//   }
// }

int main() {
  printf("=== Workspace Test ===\n");
  // KnotPointTest();
  SettingsTest();
  SolutionTest();
  DataTest();
  WorkspaceTest();
  // ProblemDataTest();
  PrintTestResult();
  return TestResult();
}
