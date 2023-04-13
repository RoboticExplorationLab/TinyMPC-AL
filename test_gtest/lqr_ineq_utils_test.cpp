// Test all inequality-related functions and tiny_RiccatiConvergence

#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

class IneqLqrTest : public testing::Test {
  public:
    const double tol = 1e-8;
    double u_max_data[NINPUTS] = {1, 1};
    double u_min_data[NINPUTS] = {-1, -1};
    double u_data[NINPUTS] = {1.1, 0.8};
    double ans[NINPUTS * 2] = {
        u_data[0] - u_max_data[0], u_data[1] - u_max_data[1],
        -u_data[0] + u_min_data[0], -u_data[1] + u_min_data[1]};
    double ineq_data[NINPUTS * 2];

    tiny_ProblemData prob;

    Matrix u;
    Matrix mat;

  private:
    void SetUp() override {
      tiny_InitProblemData(&prob);
      
      prob.nstates = NSTATES;
      prob.ninputs = NINPUTS;
      prob.nhorizon = NHORIZON;
      prob.ncstr_states = 2 * NSTATES;
      prob.ncstr_inputs = 2 * NINPUTS;
      prob.ncstr_goal = NSTATES;
      prob.u_max = slap_MatrixFromArray(NINPUTS, 1, u_max_data);
      prob.u_min = slap_MatrixFromArray(NINPUTS, 1, u_min_data);
      u = slap_MatrixFromArray(NINPUTS, 1, u_data);
      mat = slap_MatrixFromArray(NINPUTS * 2, 1, ineq_data);
    };
};

TEST_F(IneqLqrTest, IneqInputs) {
  tiny_IneqInputs(&mat, prob, u);
  // slap_PrintMatrix(mat);
  EXPECT_NEAR(SumOfSquaredError(mat.data, ans, NINPUTS * 2), 0, tol);
}

TEST_F(IneqLqrTest, IneqInputsOffset) {
  
  double ans1[NINPUTS * 2] = {u_max_data[0], u_max_data[1], -u_min_data[0],
                             -u_min_data[1]};
  tiny_IneqInputsOffset(&mat, prob);
  // slap_PrintMatrix(mat);
  EXPECT_NEAR(SumOfSquaredError(mat.data, ans1, NINPUTS * 2), 0, tol);
}

TEST_F(IneqLqrTest, IneqInputsJacobian) {
  double ans_jac[NINPUTS * 2 * NINPUTS] = {1, 0, -1, 0, 0, 1, 0, -1};
  double jac_data[NINPUTS * 2 * NINPUTS];
  Matrix jac_mat = slap_MatrixFromArray(NINPUTS * 2, NINPUTS, jac_data);
  tiny_IneqInputsJacobian(&jac_mat, prob);
  // slap_PrintMatrix(mat);
  EXPECT_NEAR(SumOfSquaredError(jac_mat.data, ans_jac, NINPUTS * 2 * NINPUTS), 0, tol);
}

TEST_F(IneqLqrTest, ActiveIneqMask) {
  double u_data1[NINPUTS] = {-2, 2};
  double ans1[NINPUTS * 2] = {u_data1[0] - u_max_data[0],    //-3
                              u_data1[1] - u_max_data[1],    // 1
                              -u_data1[0] + u_min_data[0],   // 1
                              -u_data1[1] + u_min_data[1]};  // -3
  double dual_data[NINPUTS * 2] = {1, 0, 2, 0};
  double ans2[NINPUTS * 2 * NINPUTS * 2] = {1, 0, 0, 0, 0, 1, 0, 0,
                                            0, 0, 1, 0, 0, 0, 0, 0};
  double ineq_data[NINPUTS * 2];
  double mask_data[NINPUTS * 2 * NINPUTS * 2];

  u = slap_MatrixFromArray(NINPUTS, 1, u_data1);
  Matrix ineq = slap_MatrixFromArray(NINPUTS * 2, 1, ineq_data);
  tiny_IneqInputs(&ineq, prob, u);
  Matrix input_dual = slap_MatrixFromArray(NINPUTS * 2, 1, dual_data);
  Matrix mask = slap_MatrixFromArray(NINPUTS * 2, NINPUTS * 2, mask_data);
  tiny_ActiveIneqMask(&mask, input_dual, ineq);
  // slap_PrintMatrix(mask);
  EXPECT_NEAR(SumOfSquaredError(ineq.data, ans1, NINPUTS * 2), 0, tol);
  EXPECT_NEAR(SumOfSquaredError(mask.data, ans2, NINPUTS * 2 * NINPUTS * 2), 0, tol);
}

TEST_F(IneqLqrTest, RiccatiConvergence) {
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);

  prob.nhorizon = 3;
  prob.ninputs = 2;
  double d_data[4] = {1.2, -0.3, -2.1, 3.1};
  double ans = 3.744329045369811;
  Matrix d[2];
  double* dptr = d_data;
  for (int k = 0; k < prob.nhorizon - 1; ++k) {
    d[k] = slap_MatrixFromArray(prob.ninputs, 1, dptr);
    dptr += prob.ninputs;
  }
  prob.d = d;
  double norm_d_max = tiny_RiccatiConvergence(prob);
  EXPECT_NEAR(norm_d_max, ans, 1e-6);
}

TEST_F(IneqLqrTest, ClampIneqDuals) {
  double dual_data[2] = {-2, -1};
  double new_dual_data[2] = {10, -10};
  double ans[2] = {10, 0};
  Matrix dual = slap_MatrixFromArray(2, 1, dual_data);
  Matrix new_dual = slap_MatrixFromArray(2, 1, new_dual_data);
  tiny_ClampIneqDuals(&dual, new_dual);
  EXPECT_NEAR(SumOfSquaredError(dual.data, ans, 2), 0, 1e-8);
}
