#include "simpletest.h"
#include "test_utils.h"
#include "tinympc/model_lti.h"
#include "tinympc/utils.h"
#include "bicycle_5d.h"

void tiny_SetLtiModelDims_Test() {
  // Test good 
  int NSTATES = 2;
  int NINPUTS = 1;
  tiny_LtiModel model;
  tiny_LtiModel model_not_affine;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1);

  TEST(tiny_SetLtiModelDims(&model, NSTATES, NINPUTS, 1, 0.1) == TINY_NO_ERROR);
  TEST(model.dt == (sfloat)(0.1) && model.affine == 1 && true_size == model.data_size);

  TEST(tiny_SetLtiModelDims(&model_not_affine, NSTATES, NINPUTS, 0, 0.1) == TINY_NO_ERROR);
  TEST(model_not_affine.dt == (sfloat)(0.1) && model_not_affine.affine == 0 && 
       true_size - NSTATES == model_not_affine.data_size);
  // Test bad (there are more bad cases)
  // NSTATES = 1;
  // NINPUTS = 0;
  // TEST(tiny_SetLtiModelDims(&model, NSTATES, NINPUTS) == TINY_SLAP_ERROR);
}

void tiny_InitLtiModelData_Test() {
  const sfloat tol = 1e-6;
  const int NSTATES = 2;
  const int NINPUTS = 1;
  sfloat A_data[] = {1, 0, 1, 1};  
  sfloat B_data[] = {1, 2};        
  sfloat f_data[] = {4, 5};                  
  tiny_LtiModel model;
  tiny_LtiModel model_not_affine;
  tiny_SetLtiModelDims(&model, NSTATES, NINPUTS, 1, 0.1);
  tiny_InitLtiModelData(&model, A_data, B_data, f_data);
  tiny_SetLtiModelDims(&model_not_affine, NSTATES, NINPUTS, 0, 0.1);
  tiny_InitLtiModelData(&model_not_affine, A_data, B_data, TINY_NULL);

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

  TEST(slap_IsNull(model_not_affine.f));  
}

void tiny_InitLtiModelMemory_Test() {
  const sfloat tol = 1e-6;
  const int NSTATES = 2;
  const int NINPUTS = 1;
  sfloat A_data[] = {1, 0, 1, 1};  
  sfloat B_data[] = {1, 2};        
  sfloat f_data[] = {4, 5};                  
  tiny_LtiModel model;
  tiny_SetLtiModelDims(&model, NSTATES, NINPUTS, 1, 0.1);
  sfloat data[model.data_size];
  tiny_InitLtiModelMemory(&model, data);
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

  tiny_LtiModel model_not_affine;
  tiny_SetLtiModelDims(&model_not_affine, NSTATES, NINPUTS, 0, 0.1);
  sfloat data2[model_not_affine.data_size];
  tiny_InitLtiModelMemory(&model_not_affine, data2);
  model_not_affine.A.data = (sfloat[]){1, 0, 1, 1};  // not pretty but work
  model_not_affine.B.data = B_data;
  TEST(model_not_affine.A.rows == model_not_affine.nstates);
  TEST(model_not_affine.A.cols == model_not_affine.nstates);
  TEST(SumOfSquaredError(A_data, model_not_affine.A.data, 
       model_not_affine.nstates * model_not_affine.nstates) < tol);
  TEST(model_not_affine.B.rows == model_not_affine.nstates);
  TEST(model_not_affine.B.cols == model_not_affine.ninputs);
  TEST(SumOfSquaredError(B_data, model_not_affine.B.data, 
       model_not_affine.nstates * model_not_affine.ninputs) < tol); 
  TEST(slap_IsNull(model_not_affine.f));   
}

void tiny_SetLtiModelJacFunc_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  tiny_LtiModel model;
  tiny_SetLtiModelDims(&model, NSTATES, NINPUTS, 1, 0.1);
  tiny_SetLtiModelJacFunc(&model, tiny_Bicycle5dGetJacobians);
  TEST(model.get_jacobians == tiny_Bicycle5dGetJacobians);
}

void tiny_SetLtiModelNonlFunc_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  tiny_LtiModel model;
  tiny_SetLtiModelDims(&model, NSTATES, NINPUTS, 1, 0.1);
  tiny_SetLtiModelNonlFunc(&model, tiny_Bicycle5dNonlinearDynamics);
  TEST(model.get_nonl_model == tiny_Bicycle5dNonlinearDynamics);
}

int main() {
  tiny_SetLtiModelDims_Test();
  tiny_InitLtiModelData_Test();
  tiny_InitLtiModelMemory_Test();
  tiny_SetLtiModelJacFunc_Test();
  tiny_SetLtiModelNonlFunc_Test();
  PrintTestResult();
  return TestResult();
}