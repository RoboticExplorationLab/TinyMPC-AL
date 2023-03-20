#include <gtest/gtest.h>
#include <tinympc/tinympc.h>

#include "test_utils.h"
#include "models/bicycle_5d.h"

const int NSTATES = 2;
const int NINPUTS = 1;

class ModelLtiTest : public testing::Test {
  protected:
    void SetUp() override {

    }
    
  public:
    const sfloat tol = 1e-8;
    sfloat A_data[NSTATES*NSTATES] = {1, 0, 1, 1};  
    sfloat B_data[NSTATES*NINPUTS] = {1, 2};        
    sfloat f_data[NSTATES] = {4, 5};  
};

TEST_F(ModelLtiTest, SetModelDims_Lti_Test) {
  // Test good 
  tiny_LtiModel model;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1);
  EXPECT_TRUE(tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS) == TINY_NO_ERROR 
              && true_size == model.data_size);
  // Test bad (there are more bad cases)
  // EXPECT_TRUE(tiny_SetModelDims_Lti(&model, 0, NINPUTS) == TINY_SLAP_ERROR);
}

TEST_F(ModelLtiTest, InitModelData_Lti_Test) {               
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  tiny_InitModelData_Lti(&model, A_data, B_data, f_data);

  EXPECT_TRUE(model.A.rows == model.nstates);
  EXPECT_TRUE(model.A.cols == model.nstates);
  EXPECT_TRUE(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates) < tol);
  EXPECT_TRUE(model.B.rows == model.nstates);
  EXPECT_TRUE(model.B.cols == model.ninputs);
  EXPECT_TRUE(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs) < tol);
  EXPECT_TRUE(model.f.rows == model.nstates);
  EXPECT_TRUE(model.f.cols == 1);
  EXPECT_TRUE(SumOfSquaredError(f_data, model.f.data, model.nstates) < tol);
}

TEST_F(ModelLtiTest, InitModelMemory_Lti_Test) {            
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  sfloat data[model.data_size];
  tiny_InitModelMemory_Lti(&model, data);
  model.A.data = A_data;  
  model.B.data = B_data;
  model.f.data = f_data;

  EXPECT_TRUE(model.A.rows == model.nstates);
  EXPECT_TRUE(model.A.cols == model.nstates);
  EXPECT_TRUE(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates) < tol);
  EXPECT_TRUE(model.B.rows == model.nstates);
  EXPECT_TRUE(model.B.cols == model.ninputs);
  EXPECT_TRUE(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs) < tol);
  EXPECT_TRUE(model.f.rows == model.nstates);
  EXPECT_TRUE(model.f.cols == 1);
  EXPECT_TRUE(SumOfSquaredError(f_data, model.f.data, model.nstates) < tol);
}

TEST_F(ModelLtiTest, SetModelJacFunc_Lti_Test) {
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  tiny_SetModelJacFunc_Lti(&model, tiny_Bicycle5dGetJacobians);
  EXPECT_TRUE(model.get_jacobians == tiny_Bicycle5dGetJacobians);
}

TEST_F(ModelLtiTest, SetModelNonlinear_Lti_Test) {
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  tiny_SetModelNonlinear_Lti(&model, tiny_Bicycle5dNonlinearDynamics);
  EXPECT_TRUE(model.get_nonlinear_dynamics == tiny_Bicycle5dNonlinearDynamics);
}
