#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"
#include "models/bicycle_5d.h"

const int NSTATES = 2;
const int NINPUTS = 1;
const int NHORIZON = 3;

class ModelLtvTest : public testing::Test {
  protected:
    void SetUp() override {
      tiny_MatricesfromArray(A, NSTATES, NSTATES, NHORIZON-1, A_data);
      tiny_MatricesfromArray(B, NSTATES, NINPUTS, NHORIZON-1, B_data);
      tiny_MatricesfromArray(f, NSTATES, 1, NHORIZON-1, f_data);
    }
    
  public:
    const sfloat tol = 1e-8;
    sfloat A_data[NSTATES*NSTATES*(NHORIZON-1)] = {1, 0, 1, 1, 1, 2, 3, 4};  
    sfloat B_data[NSTATES*NINPUTS*(NHORIZON-1)] = {1, 2, 5, 6};        
    sfloat f_data[NSTATES*(NHORIZON-1)] = {4, 5, 1, 1};  
    Matrix A[NHORIZON-1];
    Matrix B[NHORIZON-1];
    Matrix f[NHORIZON-1];
};

TEST_F(ModelLtvTest, SetModelDims_Ltv_Test) {
  // Test good 
  tiny_LtvModel model;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1)*(NHORIZON-1);
  EXPECT_TRUE(tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON) == TINY_NO_ERROR && true_size == model.data_size);
}

TEST_F(ModelLtvTest, InitModelData_Ltv_Test) {               
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);  
  tiny_InitModelData_Ltv(&model, A, B, f);

  for (int k = 0; k < NHORIZON-1; ++k) {
    EXPECT_TRUE(model.A[k].rows == model.nstates);
    EXPECT_TRUE(model.A[k].cols == model.nstates);
    
    EXPECT_TRUE(model.B[k].rows == model.nstates);
    EXPECT_TRUE(model.B[k].cols == model.ninputs);
    
    EXPECT_TRUE(model.f[k].rows == model.nstates);
    EXPECT_TRUE(model.f[k].cols == 1);
  }

  EXPECT_TRUE(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON - 1) < tol);
  EXPECT_TRUE(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON - 1) < tol);
  EXPECT_TRUE(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON - 1) < tol);
}

TEST_F(ModelLtvTest, InitModelMemory_Ltv_Test) {     
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  sfloat data[model.data_size];
  Matrix mats[3*(NHORIZON-1)];
  tiny_InitModelMemory_Ltv(&model, mats, data);
  model.A = A;
  model.B = B;
  model.f = f;

  for (int k = 0; k < NHORIZON-1; ++k) {
    EXPECT_TRUE(model.A[k].rows == model.nstates);
    EXPECT_TRUE(model.A[k].cols == model.nstates);
    
    EXPECT_TRUE(model.B[k].rows == model.nstates);
    EXPECT_TRUE(model.B[k].cols == model.ninputs);
    
    EXPECT_TRUE(model.f[k].rows == model.nstates);
    EXPECT_TRUE(model.f[k].cols == 1);
  }

  EXPECT_TRUE(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON - 1) < tol);
  EXPECT_TRUE(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON - 1) < tol);
  EXPECT_TRUE(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON - 1) < tol);
}

TEST_F(ModelLtvTest, SetModelJacFunc_Ltv_Test) {
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  tiny_SetModelJacFunc_Ltv(&model, tiny_Bicycle5dGetJacobians);
  EXPECT_TRUE(model.get_jacobians == tiny_Bicycle5dGetJacobians);
}

TEST_F(ModelLtvTest, SetModelNonlinear_Ltv_Test) {
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  tiny_SetModelNonlinear_Ltv(&model, tiny_Bicycle5dNonlinearDynamics);
  EXPECT_TRUE(model.get_nonlinear_dynamics == tiny_Bicycle5dNonlinearDynamics);
}
