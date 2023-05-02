#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data/forward_pass_data.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/lqr.h"
#include "tinympc/auxil.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

sfloat A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                    0.1, 0, 1, 0, 0, 0.1, 0, 1};
sfloat B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
sfloat f_data[NSTATES] = {0, 0, 0, 0};
// sfloat x0_data[NSTATES] = {5,7,2,-1.4};
sfloat Q_data[NSTATES*NSTATES] = {0};
sfloat R_data[NINPUTS*NINPUTS] = {0};
sfloat Qf_data[NSTATES] = {0};
sfloat q_data[NSTATES*(NHORIZON-1)] = {0};
sfloat r_data[NINPUTS*(NHORIZON-1)] = {0};
sfloat qf_data[NSTATES] = {0};

sfloat X_ref_data[NSTATES] = {0};
sfloat U_ref_data[NINPUTS] = {0};

void ForwardPassTest() {
  const sfloat tol = 1e-6;
  Matrix A;
  Matrix B;
  Matrix f;

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

  tiny_InitTempData(&work, temp_data);
  tiny_InitModelDataArray(&model, &A, &B, &f, A_data, B_data, f_data);

  Matrix X[NHORIZON];
  Matrix Xsln[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix X_ref[NHORIZON];
  Matrix U_ref[NHORIZON-1];
  Matrix q[NHORIZON-1];
  Matrix r[NHORIZON-1];

  sfloat* xptr = x_data;
  sfloat* xsol_ptr = xsol_data;
  sfloat* uptr = u_data;
  sfloat* Kptr = K_data;
  sfloat* dptr = d_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
      uptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;
    Xsln[i] = slap_MatrixFromArray(NSTATES, 1, xsol_ptr);
    xsol_ptr += NSTATES;
  }

  data.x0 = X[0];  // check if possible
  soln.K = K;
  soln.d = d;
  soln.U = U;
  soln.X = X;
  data.X_ref = X_ref;
  data.U_ref = U_ref;
  data.q = q;
  data.r = r;

  data.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(data.Q, 1);
  data.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(data.R, 1);
  data.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(data.Qf, 1);
  data.q[0] = slap_MatrixFromArray(NSTATES, 1, q_data);
  data.q[1] = slap_MatrixFromArray(NSTATES, 1, &q_data[NSTATES]);
  data.r[0] = slap_MatrixFromArray(NINPUTS, 1, r_data);
  data.r[1] = slap_MatrixFromArray(NINPUTS, 1, &r_data[NINPUTS]);
  data.qf = slap_MatrixFromArray(NSTATES, 1, qf_data);    
  X_ref[0] = slap_MatrixFromArray(NSTATES, 1, X_ref_data);
  X_ref[1] = slap_MatrixFromArray(NSTATES, 1, X_ref_data);
  X_ref[2] = slap_MatrixFromArray(NSTATES, 1, X_ref_data);
  U_ref[0] = slap_MatrixFromArray(NINPUTS, 1, U_ref_data);
  U_ref[1] = slap_MatrixFromArray(NINPUTS, 1, U_ref_data);
  tiny_UpdateLinearCost(&work);

  uptr = u_data;
  xsol_ptr = xsol_data;
  xptr = x_data;
  Kptr = K_data;
  dptr = d_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      TEST(U[i].rows == NINPUTS);
      TEST(U[i].cols == 1);
      TEST(SumOfSquaredError(U[i].data, uptr, NINPUTS) < tol);
      uptr += NINPUTS;
      TEST(soln.K[i].rows == NINPUTS);
      TEST(soln.K[i].cols == NSTATES);
      TEST(SumOfSquaredError(soln.K[i].data, Kptr, NINPUTS * NSTATES) < tol);
      Kptr += NINPUTS * NSTATES;
      TEST(soln.d[i].rows == NINPUTS);
      TEST(soln.d[i].cols == 1);
      TEST(SumOfSquaredError(soln.d[i].data, dptr, NINPUTS) < tol);
      dptr += NINPUTS;
    }
    TEST(X[i].rows == NSTATES);
    TEST(X[i].cols == 1);
    TEST(SumOfSquaredError(X[i].data, xptr, NSTATES) < tol);
    xptr += NSTATES;
    TEST(Xsln[i].rows == NSTATES);
    TEST(Xsln[i].cols == 1);
    TEST(SumOfSquaredError(Xsln[i].data, xsol_ptr, NSTATES) < tol);
    xsol_ptr += NSTATES;
  }

  // Include discrete dynamics test
  tiny_ForwardPass(&work);
  for (int i = 0; i < NHORIZON; ++i) {
    TEST(SumOfSquaredError(soln.X[i].data, Xsln[i].data, NSTATES) < tol);
  }
  //FIXME: cost is not exact!!
  // printf("%f\n", info.obj_pri);
  // PrintMatrix(soln.U[0]);
  // PrintMatrix(soln.U[1]);
  // PrintMatrix(soln.X[0]);PrintMatrix(soln.X[1]);PrintMatrix(soln.X[2]);
}

int main() {
  printf("=== Forward Pass Test ===\n");
  ForwardPassTest();
  PrintTestResult();
  return TestResult();
}
