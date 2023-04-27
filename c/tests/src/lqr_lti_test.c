// Test LQR
// Scenerio: Drive double integrator to arbitrary goal state.
// This one is absolutely correct

#include "tinympc/lqr_lti.h"

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/data_struct.h"
#include "tinympc/utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 100
// U, X, Psln
void LqrLtiTest() {
  sfloat A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                      0.1, 0, 1, 0, 0, 0.1, 0, 1};
  sfloat B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
  sfloat f_data[NSTATES] = {0};
  sfloat x0_data[NSTATES] = {5, 7, 2, -1.4};
  sfloat xg_data[NSTATES] = {2, 5, 3, -1};
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

  tiny_LtiModel model;
  tiny_InitLtiModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  Matrix xg = slap_MatrixFromArray(NSTATES, 1, xg_data);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];

  sfloat* Xptr = X_data;
  sfloat* Xref_ptr = Xref_data;
  sfloat* Uptr = U_data;
  sfloat* Uref_ptr = Uref_data;
  sfloat* Kptr = K_data;
  sfloat* dptr = d_data;
  sfloat* Pptr = P_data;
  sfloat* pptr = p_data;
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
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    slap_Copy(Xref[i], xg);
    Xref_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
  }
  slap_Copy(X[0], model.x0);
  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;

  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1.0);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1.0);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 1000.0);

  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  solver.reg = 1e-8;
  solver.max_outer_iters = 10;

  sfloat Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  tiny_BackwardPassLti(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLti(X, U, prob, model);

  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("\nk = \n", k);
    // tiny_Print(p[k]);
    // tiny_PrintT(Xref[k]);
    // tiny_PrintT(U[k]);
    // tiny_PrintT(X[k]);
  }
  tiny_Print(X[NHORIZON-1]);
  TEST(SumOfSquaredError(X[NHORIZON - 1].data, xg_data, NSTATES) < 1e-1);
}

int main() {
  printf("=== LQR LTI Test ===\n");
  LqrLtiTest();
  PrintTestResult();
  return TestResult();
}
