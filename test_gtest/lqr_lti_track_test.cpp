// Test tracking LQR
// Scenerio: Drive double integrator to track reference.

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"
#include "data/lqr_lti_track_data.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51

class LqrLtiTrackTest : public testing::Test {
  public:
    double A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                        0.1, 0, 1, 0, 0, 0.1, 0, 1};
    double B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
    double f_data[NSTATES] = {0};
    double x0_data[NSTATES] = {2, 6, 3, -1.5};
    double X_data[NSTATES * NHORIZON] = {0};
    double U_data[NINPUTS * (NHORIZON - 1)] = {0};
    double K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
    double d_data[NINPUTS * (NHORIZON - 1)] = {0};
    double P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
    double p_data[NSTATES * NHORIZON] = {0};
    double Q_data[NSTATES * NSTATES] = {0};
    double R_data[NINPUTS * NINPUTS] = {0};
    double Qf_data[NSTATES * NSTATES] = {0};

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

TEST_F(LqrLtiTrackTest, LqrLti) {
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
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 1000 * 1e-1);
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

  double G_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix G_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       G_temp_data);

  for (int i = 0; i < NHORIZON - 1; ++i) {
    slap_Copy(model.f, Xref[i + 1]);
    slap_MatMulAdd(model.f, model.A, Xref[i], -1, 1);
    slap_MatMulAdd(model.f, model.B, Uref[i], -1, 1);
    // tiny_Print(model.f);  // Check if reference is feasible
    // tiny_Print(slap_Transpose(Xref[i]));
  }
  tiny_BackwardPassLti(&prob, solver, model, &G_temp);
  tiny_ForwardPassLti(X, U, prob, model);
  // // tiny_AugmentedLagrangianLqr(X, U, prob, model, solver, 1);
  for (int k = 0; k < NHORIZON; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(X[k], Xref[k]));
  }
  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    EXPECT_NEAR(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES), 0, 1e-1);
  }
}
