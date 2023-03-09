#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "augmented_lagrangian_lqr.h"
#include "simpletest/simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "data/back_pass_data.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3
//U, X, Psln
void BackPassTest() {
  double A_data[NSTATES*NSTATES] = {1,0,0,0, 0,1,0,0, 0.1,0,1,0, 0,0.1,0,1};
  double B_data[NSTATES*NINPUTS] = {0.005,0,0.1,0, 0,0.005,0,0.1};
  double f_data[NSTATES] = {0};
  double x0_data[NSTATES] = {5,7,2,-1.4};
  double Xref_data[NSTATES*NHORIZON] = {0};
  double Uref_data[NINPUTS*(NHORIZON-1)] = {0};
  // double X_data[NSTATES*NHORIZON] = {0};
  // double U_data[NINPUTS*(NHORIZON-1)] = {0};
  double K_data[NINPUTS*NSTATES*(NHORIZON-1)] = {0};
  double d_data[NINPUTS*(NHORIZON-1)] = {0};
  double P_data[NSTATES*NSTATES*(NHORIZON)] = {0};
  double p_data[NSTATES*NHORIZON] = {0};
  double Q_data[NSTATES*NSTATES] = {0};
  double R_data[NINPUTS*NINPUTS] = {0};
  double Qf_data[NSTATES*NSTATES] = {0};
  double umin_data[NINPUTS] = {-2, -2};
  double umax_data[NINPUTS] = {2, 2};
  double udual_data[NINPUTS*(NHORIZON-1)] = {0};
  const double tol = 1e-8;

  tiny_LinearDiscreteModel model = kDefaultLinearDiscreteModel;
  tiny_ProblemData prob = kDefaultProblemData;
  tiny_Solver solver = kDefaultSolver;

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON-1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON-1];
  Matrix uduals[NHORIZON-1];
  Matrix K[NHORIZON-1];
  Matrix d[NHORIZON-1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];

  double* Xptr = X_data;
  double* Xref_ptr = Xref_data;
  double* Uptr = U_data;
  double* Uref_ptr = Uref_data;
  double* udual_ptr = udual_data;
  double* Kptr = K_data;
  double* dptr = d_data;
  double* Pptr = P_data;
  double* pptr = p_data;

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_MatrixFromArray(NINPUTS, 1, Uptr);
      Uptr += NINPUTS;
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
      Uref_ptr += NINPUTS;
      uduals[i] = slap_MatrixFromArray(NINPUTS, 1, udual_ptr);
      udual_ptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS*NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
    Xptr += NSTATES;    
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
    Xref_ptr += NSTATES;   
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES*NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
    pptr += NSTATES;
  }  

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.ncstr_inputs = 2*NINPUTS*(NHORIZON-1);
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 1e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);  
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 100*1e-1);
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, umax_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, umin_data);
  prob.X_ref = Xref;
  prob.U_ref = Uref;
  prob.x0 = model.x0;  
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;
  
  solver.regu = 1e-8;
  solver.input_duals = uduals;
  solver.penalty_mul = 10;

  double G_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix G_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1, 
                                      G_temp_data);
  tiny_BackwardPassLti(prob, model, solver, X, U, G_temp);
  TEST(SumOfSquaredError(d_data, dsln_data, (NHORIZON-1)*NINPUTS) < tol);
}

int main() {
  BackPassTest();
  PrintTestResult();
  return TestResult();
}
