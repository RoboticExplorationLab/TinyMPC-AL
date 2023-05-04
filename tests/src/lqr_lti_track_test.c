// Test tracking LQR
// Scenerio: Drive sfloat integrator to track reference.

#include "data/lqr_lti_track_data.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/lqr.h"
#include "tinympc/auxil.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51
// U, X, Psln
void LqrLtiTest() {
  sfloat A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                      0.1, 0, 1, 0, 0, 0.1, 0, 1};
  sfloat B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
  sfloat f_data[NSTATES] = {0};
  sfloat x0_data[NSTATES] = {2, 6, 3, -1.5};
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

  Matrix A;
  Matrix B;
  Matrix f;
  tiny_InitModelFromArray(&model, &A, &B, &f, A_data, B_data, f_data);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix q[NHORIZON-1];
  Matrix r[NHORIZON-1];

  sfloat* Xptr = X_data;
  sfloat* Xref_ptr = Xref_data;
  sfloat* Uptr = U_data;
  sfloat* Uref_ptr = Uref_data;
  sfloat* Kptr = K_data;
  sfloat* dptr = d_data;
  sfloat* Pptr = P_data;
  sfloat* pptr = p_data;
  sfloat* qptr = q_data;
  sfloat* rptr = r_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_MatrixFromArray(NINPUTS, 1, Uptr);
      // slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
      q[i] = slap_MatrixFromArray(NSTATES, 1, qptr);
      qptr += NSTATES;
      r[i] = slap_MatrixFromArray(NINPUTS, 1, rptr);
      rptr += NINPUTS;      
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    Xref_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
  }

  data.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);  // check if possible
  soln.X = X;
  soln.U = U;
  soln.K = K;
  soln.d = d;
  soln.P = P;
  soln.p = p;
  soln.U = U;
  soln.X = X;
  data.X_ref = Xref;
  data.U_ref = Uref;
  data.q = q;
  data.r = r;
  data.qf = slap_MatrixFromArray(NSTATES, 1, qf_data);  

  data.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(data.Q, 100.0);
  data.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(data.R, 1.0);
  data.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(data.Qf, 1000.0);

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

  tiny_SolveLqr(&work);

  if (0) {
    for (int k = 0; k < NHORIZON - 1; ++k) {
      printf("\n=>k = %d\n", k);
      printf("ex = %.4f\n", slap_NormedDifference(Xref[k], X[k]));
      // PrintMatrix(p[k]);
      // PrintMatrixT(Xref[k]);
      // PrintMatrixT(U[k]);
      // PrintMatrixT(X[k]);
    }
    PrintMatrixT(X[NHORIZON - 1]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    TEST(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES) < 1e-1);
  }
}

int main() {
  printf("=== LQR LTI Tracking Test ===\n");
  LqrLtiTest();
  PrintTestResult();
  return TestResult();
}
