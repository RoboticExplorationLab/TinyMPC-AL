// Test LQR 
// Scenerio: Drive double integrator to arbitrary goal state.

#include "unconstrained_lqr.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tiny_utils.h"
#include "bicycle.h"

#define H 0.1
#define NSTATES 5
#define NINPUTS 2
#define NHORIZON 101
// U, X, Psln
void LqrLtiTest() {
  double A_data[NSTATES * NSTATES * (NHORIZON-1)] = {0};
  double B_data[NSTATES * NINPUTS * (NHORIZON-1)] = {0};
  double f_data[NSTATES * (NHORIZON-1)] = {0};
  double x0_data[NSTATES] = {5, 7, 2, -1.4};
  double X_data[NSTATES * NHORIZON] = {0};
  double U_data[NINPUTS * (NHORIZON - 1)] = {0};
  double K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  double d_data[NINPUTS * (NHORIZON - 1)] = {0};
  double P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  double p_data[NSTATES * NHORIZON] = {0};
  double Q_data[NSTATES * NSTATES] = {0};
  double R_data[NINPUTS * NINPUTS] = {0};
  double Qf_data[NSTATES * NSTATES] = {0};
  double umin_data[NINPUTS] = {-2, -2};
  double umax_data[NINPUTS] = {2, 2};

  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_Solver solver;
  tiny_InitSolver(&solver);

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
  model.get_jacobians = tiny_GetJacobians;  // from Bicycle

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix A[NHORIZON-1];
  Matrix B[NHORIZON-1];
  Matrix f[NHORIZON-1];

  double* Xptr = X_data;
  double* Xref_ptr = Xref_data;
  double* Uptr = U_data;
  double* Uref_ptr = Uref_data;
  double* Kptr = K_data;
  double* dptr = d_data;
  double* Pptr = P_data;
  double* pptr = p_data;
  double* Aptr = A_data;
  double* Bptr = B_data;
  double* fptr = f_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      model.A[i] = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
      Aptr += NSTATES*NSTATES;
      model.B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
      Bptr += NSTATES*NINPUTS;
      model.f[i] = slap_MatrixFromArray(NSTATES, 1, f_data);
      fptr += NSTATES;
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
    Xref_ptr += NSTATES;
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
  }
 
  slap_MatrixCopy(X[0], model.x0);
  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 2 * NINPUTS * (NHORIZON - 1);
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 1000 * 1e-1);
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, umax_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, umin_data);
  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  double Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  for (int i = 0; i < NHORIZON-1; ++i) {  
    model.get_jacobians(&(model.A[i]), &(model.B[i]), prob.X_ref[i], prob.U_ref[i]);
  }                                     

  tiny_BackwardPassLtv(&prob, solver, model, Q_temp);
  tiny_ForwardPassLtv(X, U, prob, model);

  // tiny_AugmentedLagrangianLqr(X, U, prob, model, solver, 1);
  for (int k = 0; k < NHORIZON-1; ++k) {
    tiny_Print(U[k]);
  }
  tiny_Print(X[NHORIZON-1]);
  TEST(SumOfSquaredError(X[NHORIZON - 1].data, Xref[NHORIZON-1].data, NSTATES) < 1e-1);
}

int main() {
  LqrLtiTest();
  PrintTestResult();
  return TestResult();
}
