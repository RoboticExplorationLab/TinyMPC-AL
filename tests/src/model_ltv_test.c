#include "simpletest.h"
#include "test_utils.h"
#include "tinympc/model_ltv.h"
#include "tinympc/utils.h"
#include "bicycle_5d.h"

void tiny_SetModelDims_Ltv_Test() {
  // Test good 
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_LtvModel model;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1)*(NHORIZON-1);
  TEST(tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON) == TINY_NO_ERROR &&
       true_size == model.data_size);
  // Test bad (there are more bad cases)
  // NSTATES = 1;
  // NINPUTS = 0;
  // int NHORIZON = 0;
  // TEST(tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON) == TINY_SLAP_ERROR);
}

void tiny_InitModelData_Ltv_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 1;
  const int NINPUTS = 1;
  const int NHORIZON = 2;
  sfloat A_data[] = {1, 0};  
  sfloat B_data[] = {1, 2};        
  sfloat f_data[] = {4, 5};                  
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  tiny_InitModelData_Ltv(&model, A_data, B_data, f_data);

  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(model.A[k].rows == model.nstates);
    TEST(model.A[k].cols == model.nstates);
    TEST(SumOfSquaredError(A_data[k], model.A[k].data, model.nstates * model.nstates) <
        tol);
    TEST(model.B[k].rows == model.nstates);
    TEST(model.B[k].cols == model.ninputs);
    TEST(SumOfSquaredError(B_data[k], model.B[k].data, model.nstates * model.ninputs) <
        tol);
    TEST(model.f[k].rows == model.nstates);
    TEST(model.f[k].cols == 1);
    TEST(SumOfSquaredError(f_data[k], model.f[k].data, model.nstates) < tol);
  }
}

// void tiny_InitModelMemory_Ltv_Test() {
//   const sfloat tol = 1e-8;
//   const int NSTATES = 2;
//   const int NINPUTS = 1;
//   sfloat A_data[] = {1, 0, 1, 1};  
//   sfloat B_data[] = {1, 2};        
//   sfloat f_data[2] = {4, 5};                  
//   tiny_LtvModel model;
//   tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS);
//   sfloat data[model.data_size];
//   tiny_InitModelMemory_Ltv(&model, data);
//   model.A.data = (sfloat[]){1, 0, 1, 1};  // not pretty but work
//   model.B.data = B_data;
//   model.f.data = (sfloat[]){4, 5};

//   TEST(model.A.rows == model.nstates);
//   TEST(model.A.cols == model.nstates);
//   TEST(SumOfSquaredError(A_data, model.A.data, model.nstates * model.nstates) <
//        tol);
//   TEST(model.B.rows == model.nstates);
//   TEST(model.B.cols == model.ninputs);
//   TEST(SumOfSquaredError(B_data, model.B.data, model.nstates * model.ninputs) <
//        tol);
//   TEST(model.f.rows == model.nstates);
//   TEST(model.f.cols == 1);
//   TEST(SumOfSquaredError(f_data, model.f.data, model.nstates) < tol);
// }

void tiny_SetModelJacFunc_Ltv_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  tiny_SetModelJacFunc_Ltv(&model, tiny_Bicycle5dGetJacobians);
  TEST(model.get_jacobians == tiny_Bicycle5dGetJacobians);
}

void tiny_SetModelNonlinear_Ltv_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  tiny_SetModelNonlinear_Ltv(&model, tiny_Bicycle5dNonlinearDynamics);
  TEST(model.get_nonlinear_dynamics == tiny_Bicycle5dNonlinearDynamics);
}

int main() {
  tiny_SetModelDims_Ltv_Test();
  // tiny_InitModelData_Ltv_Test();
  // tiny_InitModelMemory_Ltv_Test();
  tiny_SetModelJacFunc_Ltv_Test();
  tiny_SetModelNonlinear_Ltv_Test();
  PrintTestResult();
  return TestResult();
}
