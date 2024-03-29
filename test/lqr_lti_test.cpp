// Test LQR
// Scenerio: Drive double integrator to arbitrary goal state.

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51

class LqrLtiTest : public testing::Test {
  public:
    double A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                        0.1, 0, 1, 0, 0, 0.1, 0, 1};
    double B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
    double f_data[NSTATES] = {0};
    double x0_data[NSTATES] = {5, 7, 2, -1.4};
    double xg_data[NSTATES] = {2, 5, 3, -1};
    double Xref_data[NSTATES * NHORIZON] = {0};
    double Uref_data[NINPUTS * (NHORIZON - 1)] = {0};
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
    
    tiny_LtiModel model;
    tiny_ProblemData prob;
    tiny_Solver solver;
      
    Matrix X[NHORIZON];
    Matrix U[NHORIZON - 1];
    Matrix Xref[NHORIZON];
    Matrix Uref[NHORIZON - 1];
    Matrix K[NHORIZON - 1];
    Matrix d[NHORIZON - 1];
    Matrix P[NHORIZON];
    Matrix p[NHORIZON];

    double* Xptr = X_data;
    double* Xref_ptr = Xref_data;
    double* Uptr = U_data;
    double* Uref_ptr = Uref_data;
    double* Kptr = K_data;
    double* dptr = d_data;
    double* Pptr = P_data;
    double* pptr = p_data;

    Matrix Q_temp;

  private:
    void SetUp() override {
      tiny_InitLtiModel(&model);
      tiny_InitProblemData(&prob);
      tiny_InitSolver(&solver);

      model.ninputs = NSTATES;
      model.nstates = NINPUTS;
      model.A = slap_MatrixFromArray(NSTATES, NSTATES, A_data);
      model.B = slap_MatrixFromArray(NSTATES, NINPUTS, B_data);
      model.f = slap_MatrixFromArray(NSTATES, 1, f_data);
      model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);

    }
};

// U, X, Psln
TEST_F(LqrLtiTest, LqrLti) {
  Matrix xg = slap_MatrixFromArray(NSTATES, 1, xg_data);
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

  solver.reg = 1e-8;
  solver.penalty_mul = 10;
  solver.max_primal_iters = 1;
  
  double Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Q_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                      Q_temp_data);
  
  tiny_BackwardPassLti(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLti(X, U, prob, model);

  // tiny_AugmentedLagrangianLqr(X, U, prob, model, solver, 1);
  // for (int k = 0; k < NHORIZON - 1; ++k) {
  //   // tiny_Print(U[k]);
  // }
  // tiny_Print(X[NHORIZON-1]);
  EXPECT_NEAR(SumOfSquaredError(X[NHORIZON - 1].data, xg_data, NSTATES), 0, 1e-1);
}