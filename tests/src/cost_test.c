#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/cost_lqr.h"
#include "tinympc/model.h"
#include "tinympc/auxil.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 2

sfloat x_data[NSTATES] = {1.1, 1.2, 1.3, -4.3};
sfloat u_data[NINPUTS] = {-2.1, 1.1};
sfloat x_ref_data[NSTATES * NHORIZON] = {1.1, 1.2, 1.3, -4.2,
                                         1.2, 1.3, 1.3, -4.3};
sfloat u_ref_data[NINPUTS * (NHORIZON - 1)] = {-2.1, 1.4};
sfloat Q_data[NSTATES * NSTATES] = {0};  // NOLINT
sfloat R_data[NINPUTS * NINPUTS] = {0};  // NOLINT
sfloat q_data[NSTATES*NHORIZON] = {0};            // NOLINT
sfloat qf_data[NSTATES] ={0};
sfloat r_data[NINPUTS*(NHORIZON-1)] = {0};            // NOLINT
sfloat Qf_data[NSTATES * NSTATES] = {0};
sfloat ans_stage[2] = {0.04549999999999994, 0.1314999999999999};
sfloat ans_term = 0.0049999999999999975;
sfloat ans_gradx[NSTATES] = {-0.11, -0.12, -0.13, 0.42};
sfloat ans_gradu[NINPUTS] = {0.21, -0.14};
sfloat ans_gradxf[NSTATES] = {-0.6, -0.65, -0.65, 2.15};

void AddCostTest() {
  const sfloat tol = 1e-6;
  Matrix U_ref[NHORIZON-1];
  Matrix X_ref[NHORIZON];
  Matrix U[NHORIZON-1];
  Matrix X[NHORIZON];
  sfloat* uptr = u_ref_data;
  sfloat* xptr = x_ref_data;
  for (int i = 0; i < NHORIZON; ++i) {
    U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
    uptr += NINPUTS;
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
  }
  Matrix x = slap_MatrixFromArray(NSTATES, 1, x_data);
  Matrix u = slap_MatrixFromArray(NINPUTS, 1, u_data);

  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
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

  data.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(data.Q, 0.1);
  data.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(data.R, 1);
  data.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(data.Qf, 0.5);
  data.X_ref = X_ref;
  data.U_ref = U_ref;

  soln.X = X;
  soln.U = U;

  work.info->obj_pri = 0.0;
  for (int k = 0; k < NHORIZON - 1; ++k) {
    soln.X[k] = x;
    soln.U[k] = u;
    tiny_AddStageCost(&work, k);
    TESTAPPROX(work.info->obj_pri, ans_stage[k], tol);
  }
  soln.X[NHORIZON - 1] = x;
  work.info->obj_pri = 0.0;
  tiny_AddTerminalCost(&work);
  TESTAPPROX(work.info->obj_pri, ans_term, tol);
}

void ExpandCostTest() {
  const sfloat tol = 1e-6;
  Matrix U[NHORIZON-1];
  Matrix X[NHORIZON];
  Matrix U_ref[NHORIZON-1];
  Matrix X_ref[NHORIZON];
  Matrix q[NHORIZON];
  Matrix r[NHORIZON - 1];

  sfloat* uptr = u_ref_data;
  sfloat* xptr = x_ref_data;
  for (int i = 0; i < NHORIZON; ++i) {
    U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
    uptr += NINPUTS;
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
  }

  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
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
  data.q = q;
  data.r = r;
  data.qf = slap_MatrixFromArray(NSTATES, 1, ans_gradxf);
  soln.X = X;
  soln.U = U;

  data.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(data.Q, 0.1);
  data.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(data.R, 0.1);
  data.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(data.Qf, 0.5);
  data.q[0] = slap_MatrixFromArray(NSTATES, 1, q_data);
  data.q[1] = slap_MatrixFromArray(NSTATES, 1, &q_data[NSTATES]);
  data.r[0] = slap_MatrixFromArray(NINPUTS, 1, r_data);
  data.qf = slap_MatrixFromArray(NSTATES, 1, qf_data);    

  data.X_ref = X_ref;
  data.U_ref = U_ref;
  tiny_UpdateLinearCost(&work);
  
  tiny_ExpandStageCost(&work, 0);
  TEST(SumOfSquaredError(work.Qx.data, ans_gradx, NINPUTS * NINPUTS) < tol);
  TEST(SumOfSquaredError(work.Qu.data, ans_gradu, NINPUTS) < tol);
  // FIXME: problem with sub-matrices
  // tiny_ExpandTerminalCost(&work);  // have to allocate P, p first
  // TEST(SumOfSquaredError(work.Qxx.data, data.Qf.data, NSTATES * NSTATES) < tol);
  // TEST(SumOfSquaredError(work.soln->p[NHORIZON-1].data, ans_gradxf, NSTATES) < tol);
}

int main() {
  printf("=== Cost Test ===\n");
  AddCostTest();
  ExpandCostTest();
  PrintTestResult();
  return TestResult();
}
