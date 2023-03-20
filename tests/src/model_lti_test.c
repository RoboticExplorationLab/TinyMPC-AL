#include "simpletest.h"
#include "test_utils.h"
#include "tinympc/model_lti.h"
#include "tinympc/utils.h"
#include "bicycle_5d.h"

void tiny_SetModelDims_Lti_Test() {
  // Test good 
  int NSTATES = 2;
  int NINPUTS = 1;
  tiny_LtiModel model;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1);
  TEST(tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS) == TINY_NO_ERROR &&
       true_size == model.data_size);
  // Test bad (there are more bad cases)
  // NSTATES = 1;
  // NINPUTS = 0;
  // TEST(tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS) == TINY_SLAP_ERROR);
}

void tiny_InitModelData_Lti_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 2;
  const int NINPUTS = 1;
  sfloat A_data[] = {1, 0, 1, 1};  
  sfloat B_data[] = {1, 2};        
  sfloat f_data[] = {4, 5};                  
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  tiny_InitModelData_Lti(&model, A_data, B_data, f_data);

  TEST(model.A.rows == model.nstates);
  TEST(model.A.cols == model.nstates);
  TEST(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates) <
       tol);
  TEST(model.B.rows == model.nstates);
  TEST(model.B.cols == model.ninputs);
  TEST(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs) <
       tol);
  TEST(model.f.rows == model.nstates);
  TEST(model.f.cols == 1);
  TEST(SumOfSquaredError(f_data, model.f.data, model.nstates) < tol);
}

void tiny_InitModelMemory_Lti_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 2;
  const int NINPUTS = 1;
  sfloat A_data[] = {1, 0, 1, 1};  
  sfloat B_data[] = {1, 2};        
  sfloat f_data[] = {4, 5};                  
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  sfloat data[model.data_size];
  tiny_InitModelMemory_Lti(&model, data);
  model.A.data = (sfloat[]){1, 0, 1, 1};  // not pretty but work
  model.B.data = B_data;
  model.f.data = (sfloat[]){4, 5};

  TEST(model.A.rows == model.nstates);
  TEST(model.A.cols == model.nstates);
  TEST(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates) <
       tol);
  TEST(model.B.rows == model.nstates);
  TEST(model.B.cols == model.ninputs);
  TEST(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs) <
       tol);
  TEST(model.f.rows == model.nstates);
  TEST(model.f.cols == 1);
  TEST(SumOfSquaredError(f_data, model.f.data, model.nstates) < tol);
}

void tiny_SetModelJacFunc_Lti_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  tiny_SetModelJacFunc_Lti(&model, tiny_Bicycle5dGetJacobians);
  TEST(model.get_jacobians == tiny_Bicycle5dGetJacobians);
}

void tiny_SetModelNonlinear_Lti_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  tiny_LtiModel model;
  tiny_SetModelDims_Lti(&model, NSTATES, NINPUTS);
  tiny_SetModelNonlinear_Lti(&model, tiny_Bicycle5dNonlinearDynamics);
  TEST(model.get_nonlinear_dynamics == tiny_Bicycle5dNonlinearDynamics);
}

int main() {
  tiny_SetModelDims_Lti_Test();
  tiny_InitModelData_Lti_Test();
  tiny_InitModelMemory_Lti_Test();
  tiny_SetModelJacFunc_Lti_Test();
  tiny_SetModelNonlinear_Lti_Test();
  PrintTestResult();
  return TestResult();
}
