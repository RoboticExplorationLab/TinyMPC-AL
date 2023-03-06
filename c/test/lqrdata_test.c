#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "augmented_lagrangian_lqr.h"
#include "simpletest/simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"

#define NSTATES 2
#define NINPUTS 1
#define NHORIZON 3

// int nx = 2;
// int nu = 1;
// int nh = 2;
double dt = 0.1;
double A_data[NSTATES*NSTATES] = {1, 0, 1, 1};  // NOLINT
double B_data[NSTATES*NINPUTS] = {1, 2};        // NOLINT
double f_data[NSTATES] = {4, 5};        // NOLINT
double t = 1.0;
double Q_data[NSTATES*NSTATES] = {1, 0, 0, 1};  // NOLINT
double R_data[NINPUTS*NINPUTS] = {1};         // NOLINT
double q_data[NSTATES] = {0.1, 0.2};    // NOLINT
double r_data[NINPUTS] = {-0.6};        // NOLINT
double Qf_data[NSTATES*NSTATES] = {0};
double u_max_data[NINPUTS] = {1.1};
double u_min_data[NINPUTS] = {-1.1};
double x_max_data[NSTATES] = {1.6, 1.7};
double x_min_data[NSTATES] = {-1.6, -1.7};
double x_ref_data[NSTATES*NHORIZON] = {0.2, 1.1,  2.5, 3.7,  2.1, 4.5};
double u_ref_data[NINPUTS*(NHORIZON-1)] = {1, 2};
double x0_data[NSTATES] = {0.1, 0.2};
double u0_data[NSTATES] = {1.6};
double Kd_data[NINPUTS*(NSTATES+1)*(NHORIZON-1)] = {0};
double Pp_data[NSTATES*(NSTATES+1)*NHORIZON] = {0};
double regu = 1e-8;
double input_duals_data[2*NINPUTS*(NHORIZON-1)] = {1, 2,  3, 4};
double state_duals_data[2*NSTATES*(NHORIZON)] = {1,2,3,4, 5,6,7,8, 2,4,5,6};
double goal_duals_data[NSTATES] = {1,2};
double penalty_min = 1;
double penalty_max = 1e8;
double penalty_mul = 10;

void LinearDiscreteModelTest() {
  const double tol = 1e-8;
  tiny_LinearDiscreteModel model = kDefaultLinearDiscreteModel;
  model.nstates = NSTATES;
  model.ninputs = NINPUTS;
  model.dt = dt;
  model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
  model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);  
  model.f = slap_MatrixFromArray(NSTATES, 1, f_data);  
  model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data); 

  TEST(model.dt == dt);
  TEST(model.A.rows == model.nstates);
  TEST(model.A.cols == model.nstates);
  TEST(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates) < tol);
  TEST(model.B.rows == model.nstates);
  TEST(model.B.cols == model.ninputs);
  TEST(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs) < tol);
  TEST(model.f.rows == model.nstates);
  TEST(model.f.cols == 1);  
  TEST(SumOfSquaredError(f_data, model.f.data, model.nstates) < tol);  
  TEST(model.x0.rows == model.nstates);
  TEST(model.x0.cols == 1);  
  TEST(SumOfSquaredError(x0_data, model.x0.data, model.nstates) < tol);    
}

void KnotPointTest() {
  const double tol = 1e-8;
  tiny_KnotPoint z = kDefaultKnotPoint;
  tiny_KnotPoint Z[NHORIZON] = {kDefaultKnotPoint};
  z.dt = dt;
  z.t = t;
  z.x = slap_MatrixFromArray(NSTATES, 1, x0_data); 
  z.u = slap_MatrixFromArray(NINPUTS, 1, u0_data); 

  TEST(z.dt == dt);
  TEST(z.t == t);
  TEST(z.x.rows == NSTATES);
  TEST(z.x.cols == 1);
  TEST(SumOfSquaredError(x0_data, z.x.data, z.x.rows) < tol);
  TEST(z.u.rows == NINPUTS);
  TEST(z.u.cols == 1);  
  TEST(SumOfSquaredError(u0_data, z.u.data, z.u.rows) < tol);

  for (int i = 0; i < NHORIZON; ++i) {
    Z[i].dt = dt;
    Z[i].t = t + i;
    Z[i].x = slap_MatrixFromArray(NSTATES, 1, x0_data); 
    Z[i].u = slap_MatrixFromArray(NINPUTS, 1, u0_data); 
  }

  for (int i = 0; i < NHORIZON; ++i) {
    TEST(Z[i].dt == dt);
    TEST(Z[i].t == (t + i));
    TEST(Z[i].x.rows == NSTATES);
    TEST(Z[i].x.cols == 1);
    TEST(SumOfSquaredError(Z[i].x.data, z.x.data, z.x.rows) < tol);
    TEST(Z[i].u.rows == NINPUTS);
    TEST(Z[i].u.cols == 1);  
    TEST(SumOfSquaredError(Z[i].u.data, z.u.data, z.u.rows) < tol);
  }
}

void SolverTest() {
  const double tol = 1e-8;
  Matrix input_duals[NHORIZON-1];
  Matrix state_duals[NHORIZON];
  double* input_dual_ptr = input_duals_data;
  double* state_dual_ptr = state_duals_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      input_duals[i] = slap_MatrixFromArray(NINPUTS, 1, input_dual_ptr);
      input_dual_ptr += NINPUTS;
    }
    state_duals[i] = slap_MatrixFromArray(NSTATES, 1, state_dual_ptr);
    state_dual_ptr += NSTATES;    
  }

  tiny_Solver solver = kDefaultSolver;
  solver.regu = regu;
  solver.input_duals = input_duals;
  solver.state_duals = state_duals;
  solver.goal_dual = slap_MatrixFromArray(NSTATES, 1, goal_duals_data);
  solver.penalty_max = penalty_max;
  solver.penalty_min = penalty_min;
  solver.penalty_mul = penalty_mul;

  TEST(solver.regu == regu);
  TEST(solver.penalty_max == penalty_max);
  TEST(solver.penalty_min == penalty_min);
  TEST(solver.penalty_mul == penalty_mul);
  TEST(SumOfSquaredError(solver.goal_dual.data, goal_duals_data, NSTATES) < tol);
  input_dual_ptr = input_duals_data;
  state_dual_ptr = state_duals_data;  
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      TEST(solver.input_duals[i].rows == NINPUTS);
      TEST(solver.input_duals[i].cols == 1);  
      TEST(SumOfSquaredError(solver.input_duals[i].data, input_dual_ptr, NINPUTS) < tol);
      input_dual_ptr += NINPUTS;
    }
      TEST(solver.state_duals[i].rows == NSTATES);
      TEST(solver.state_duals[i].cols == 1);  
      TEST(SumOfSquaredError(solver.state_duals[i].data, state_dual_ptr, NSTATES) < tol);
      state_dual_ptr += NSTATES; 
  }  
}

void ProblemDataTest() {
  const double tol = 1e-8;
  Matrix U_ref[NHORIZON-1];
  Matrix X_ref[NHORIZON];
  Matrix K[NHORIZON-1];
  Matrix d[NHORIZON-1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  double* uptr = u_ref_data;
  double* xptr = x_ref_data;  
  double* Kd_ptr = Kd_data;
  double* Pp_ptr = Pp_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U_ref[i] = slap_MatrixFromArray(NINPUTS, 1, uptr);
      uptr += NINPUTS;
      K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kd_ptr);
      Kd_ptr += NINPUTS*NSTATES;
      d[i] = slap_MatrixFromArray(NINPUTS, 1, Kd_ptr);
      Kd_ptr += NINPUTS;
    }
    X_ref[i] = slap_MatrixFromArray(NSTATES, 1, xptr);
    xptr += NSTATES;    
    P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pp_ptr);
    Pp_ptr += NSTATES * NSTATES;
    p[i] = slap_MatrixFromArray(NSTATES, 1, Pp_ptr);
    Pp_ptr += NSTATES;
  }

  tiny_ProblemData prob = kDefaultProblemData;
  prob.nstates = NSTATES;
  prob.ninputs = NINPUTS;
  prob.nhorizon = NHORIZON;
  prob.ncstr_goal = NSTATES;
  prob.ncstr_inputs = 2*NINPUTS*(NHORIZON-1);
  prob.ncstr_states = 2*NSTATES*NHORIZON;
  prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  prob.q = slap_MatrixFromArray(NSTATES, 1, q_data);
  prob.r = slap_MatrixFromArray(NINPUTS, 1, r_data);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
  prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
  prob.x_max = slap_MatrixFromArray(NSTATES, 1, x_max_data);
  prob.x_min = slap_MatrixFromArray(NSTATES, 1, x_min_data);
  prob.X_ref = X_ref;
  prob.U_ref = U_ref;
  prob.dt = dt;
  prob.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data); 
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  TEST(prob.nstates == NSTATES);
  TEST(prob.ninputs == NINPUTS);
  TEST(prob.nhorizon == NHORIZON);
  TEST(prob.ncstr_inputs == 2*NINPUTS*(NHORIZON-1));  
  TEST(prob.ncstr_states == 2*NSTATES*NHORIZON);
  TEST(prob.dt == dt);
  TEST(SumOfSquaredError(prob.Q.data, Q_data, NSTATES*NSTATES) < tol);
  TEST(SumOfSquaredError(prob.R.data, R_data, NINPUTS*NINPUTS) < tol);
  TEST(SumOfSquaredError(prob.q.data, q_data, NSTATES) < tol);
  TEST(SumOfSquaredError(prob.r.data, r_data, NINPUTS) < tol);
  TEST(SumOfSquaredError(prob.Qf.data, Qf_data, NSTATES*NSTATES) < tol);
  TEST(SumOfSquaredError(prob.u_max.data, u_max_data, NINPUTS) < tol);
  TEST(SumOfSquaredError(prob.u_min.data, u_min_data, NINPUTS) < tol);
  TEST(SumOfSquaredError(prob.x_max.data, x_max_data, NSTATES) < tol);
  TEST(SumOfSquaredError(prob.x_min.data, x_min_data, NSTATES) < tol);
  TEST(SumOfSquaredError(prob.x0.data, x0_data, NSTATES) < tol);
  uptr = u_ref_data;
  xptr = x_ref_data;  
  Kd_ptr = Kd_data;
  Pp_ptr = Pp_data;
  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      TEST(prob.U_ref[i].rows == NINPUTS);
      TEST(prob.U_ref[i].cols == 1);  
      TEST(SumOfSquaredError(prob.U_ref[i].data, uptr, NINPUTS) < tol);
      uptr += NINPUTS;
      TEST(prob.K[i].rows == NINPUTS);
      TEST(prob.K[i].cols == NSTATES);  
      TEST(SumOfSquaredError(prob.K[i].data, Kd_ptr, NINPUTS*NSTATES) < tol);
      Kd_ptr += NINPUTS*NSTATES;
      TEST(prob.d[i].rows == NINPUTS);
      TEST(prob.d[i].cols == 1);  
      TEST(SumOfSquaredError(prob.d[i].data, Kd_ptr, NINPUTS) < tol);
      Kd_ptr += NINPUTS;    
    }
    TEST(prob.X_ref[i].rows == NSTATES);
    TEST(prob.X_ref[i].cols == 1);  
    TEST(SumOfSquaredError(prob.X_ref[i].data, xptr, NSTATES) < tol);
    xptr += NSTATES;
    TEST(prob.P[i].rows == NSTATES);
    TEST(prob.P[i].cols == NSTATES);  
    TEST(SumOfSquaredError(prob.P[i].data, Pp_ptr, NSTATES*NSTATES) < tol);
    Pp_ptr += NSTATES*NSTATES;
    TEST(prob.p[i].rows == NSTATES);
    TEST(prob.p[i].cols == 1);  
    TEST(SumOfSquaredError(prob.p[i].data, Pp_ptr, NSTATES) < tol);
    Pp_ptr += NSTATES;    
  }  
}

int main() {
  LinearDiscreteModelTest();
  KnotPointTest();
  SolverTest();
  ProblemDataTest();
  PrintTestResult();
  return TestResult();
}
