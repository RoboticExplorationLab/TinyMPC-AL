// Test TVLQR
// Scenerio: Drive bicycle to track references.

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"
#include "models/bicycle_5d.h"

#define H 0.1
#define NSTATES 5
#define NINPUTS 2
#define NHORIZON 81
// NO GRADIENT VANISHING/EXPLOSION WHEN NHORIZON = 101

class LqrLtvTrackTest : public testing::Test {
  public:
    double x0_data[NSTATES] = {1, -1, 0, 0, 0};
    double xg_data[NSTATES] = {0};
    double ug_data[NINPUTS] = {0};
    double Q_data[NSTATES * NSTATES] = {0};
    double R_data[NINPUTS * NINPUTS] = {0};
    double Qf_data[NSTATES * NSTATES] = {0};

    Matrix X[NHORIZON];
    Matrix U[NHORIZON - 1];
    Matrix Xref[NHORIZON];
    Matrix Uref[NHORIZON - 1];
    Matrix Xnom[NHORIZON];
    Matrix Unom[NHORIZON - 1];
    Matrix K[NHORIZON - 1];
    Matrix d[NHORIZON - 1];
    Matrix P[NHORIZON];
    Matrix p[NHORIZON];
    Matrix A[NHORIZON - 1];
    Matrix B[NHORIZON - 1];
    Matrix f[NHORIZON - 1];
    
    double X_data[NSTATES * NHORIZON] = {0};
    double U_data[NINPUTS * (NHORIZON - 1)] = {0};
    double K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
    double d_data[NINPUTS * (NHORIZON - 1)] = {0};
    double P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
    double p_data[NSTATES * NHORIZON] = {0};
    double A_data[NSTATES * NSTATES * (NHORIZON - 1)] = {0};
    double B_data[NSTATES * NINPUTS * (NHORIZON - 1)] = {0};
    double f_data[NSTATES * (NHORIZON - 1)] = {0};
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

    tiny_LtvModel model;
    tiny_ProblemData prob;
    tiny_Solver solver;

  private:
    void SetUp() override {
      tiny_InitLtvModel(&model);
      tiny_InitProblemData(&prob);
      tiny_InitSolver(&solver);

      for (int i = 0; i < NHORIZON; ++i) {
        if (i < NHORIZON - 1) {
          A[i] = slap_MatrixFromArray(NSTATES, NSTATES, Aptr);
          Aptr += NSTATES * NSTATES;
          B[i] = slap_MatrixFromArray(NSTATES, NINPUTS, Bptr);
          Bptr += NSTATES * NINPUTS;
          f[i] = slap_MatrixFromArray(NSTATES, 1, fptr);
          fptr += NSTATES;
          U[i] = slap_MatrixFromArray(NINPUTS, 1, Uptr);
          // slap_SetConst(U[i], 0.01);
          Uptr += NINPUTS;
          Unom[i] = slap_MatrixFromArray(NINPUTS, 1, Uref_ptr);
          Uref_ptr += NINPUTS;
          Uref[i] = slap_MatrixFromArray(NINPUTS, 1, ug_data);
          K[i] = slap_MatrixFromArray(NINPUTS, NSTATES, Kptr);
          Kptr += NINPUTS * NSTATES;
          d[i] = slap_MatrixFromArray(NINPUTS, 1, dptr);
          dptr += NINPUTS;
        }
        X[i] = slap_MatrixFromArray(NSTATES, 1, Xptr);
        Xptr += NSTATES;
        Xnom[i] = slap_MatrixFromArray(NSTATES, 1, Xref_ptr);
        Xref_ptr += NSTATES;
        Xref[i] = slap_MatrixFromArray(NSTATES, 1, xg_data);
        P[i] = slap_MatrixFromArray(NSTATES, NSTATES, Pptr);
        Pptr += NSTATES * NSTATES;
        p[i] = slap_MatrixFromArray(NSTATES, 1, pptr);
        pptr += NSTATES;
      }

      model.ninputs = NSTATES;
      model.nstates = NINPUTS;
      model.x0 = slap_MatrixFromArray(NSTATES, 1, x0_data);
      model.get_jacobians = tiny_Bicycle5dGetJacobians;  // from Bicycle
      model.A = A;
      model.B = B;
      model.f = f;
      slap_Copy(X[0], model.x0);

      prob.ninputs = NINPUTS;
      prob.nstates = NSTATES;
      prob.nhorizon = NHORIZON;
      prob.Q = slap_MatrixFromArray(NSTATES, NSTATES, Q_data);
      slap_SetIdentity(prob.Q, 10e-1);
      prob.R = slap_MatrixFromArray(NINPUTS, NINPUTS, R_data);
      slap_SetIdentity(prob.R, 1e-1);
      prob.Qf = slap_MatrixFromArray(NSTATES, NSTATES, Qf_data);
      slap_SetIdentity(prob.Qf, 10e-1);
      prob.X_ref = Xref;
      prob.U_ref = Uref;
      prob.x0 = model.x0;
      prob.K = K;
      prob.d = d;
      prob.P = P;
      prob.p = p;
    }
}


TEST_F(LqrLtvTrackTest, DeltaLqrLtv) {
  double Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  // DELTA FORMULATION: Xref in objective and Xnom linearization
  // Compute and store A, B offline
  for (int i = 0; i < NHORIZON - 1; ++i) {
    model.get_jacobians(&(model.A[i]), &(model.B[i]), Xnom[i], Unom[i]);
    tiny_Bicycle5dNonlinearDynamics(&(model.f[i]), Xnom[i], Unom[i]);
    slap_MatrixAddition(model.f[i], model.f[i], Xnom[i + 1], -1);  // = 0
    // printf("%f\n",slap_NormTwoSquared(model.f[i]));
  }

  tiny_BackwardPassLtv(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLtv(X, U, prob, model);  // delta_x and delta_u

  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(X[k], Xref[k]));
    // tiny_NonlinearDynamics(&X[k+1], X[k], Uref[k]);
    // tiny_Print(Uref[k]);
    // tiny_Print(model.B[k]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    EXPECT_LT(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES), 0.5);
  }
}

TEST_F(LqrLtvTrackTest, AbsLqrLtv) {
  double Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_MatrixFromArray(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  // Absolute formulation
  // Compute and store A, B offline
  for (int i = 0; i < NHORIZON - 1; ++i) {
    model.get_jacobians(&(model.A[i]), &(model.B[i]), prob.X_ref[i],
                        prob.U_ref[i]);
    tiny_Bicycle5dNonlinearDynamics(&(model.f[i]), Xref[i], Uref[i]);
    slap_MatMulAdd(model.f[i], model.A[i], Xref[i], -1, 1);
    slap_MatMulAdd(model.f[i], model.B[i], Uref[i], -1, 1);
    // printf("%f\n",slap_NormTwoSquared(model.f[i]));
  }

  tiny_BackwardPassLtv(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLtv(X, U, prob, model);

  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_MatrixNormedDifference(X[k], Xref[k]));
    // tiny_NonlinearDynamics(&X[k+1], X[k], Uref[k]);
    // tiny_Print(slap_Transpose(X[k]));
    // tiny_Print(model.B[k]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    EXPECT_LT(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES), 0.5);
  }
}