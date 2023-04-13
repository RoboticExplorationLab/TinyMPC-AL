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

void tiny_InitModelDataArray_Ltv_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 2;
  const int NINPUTS = 2;
  const int NHORIZON = 3;
  sfloat A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
  sfloat B_data[] = {1.1, 2.2, 1.3, 2.3};        
  sfloat f_data[] = {2.1, 1.2, 3.3, 1.3};     
  Matrix A[NHORIZON-1];
  Matrix B[NHORIZON-1];     
  Matrix f[NHORIZON-1];        
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  tiny_InitModelDataArray_Ltv(&model, A, B, f, A_data, B_data, f_data);
  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(model.A[k].rows == model.nstates);
    TEST(model.A[k].cols == model.nstates);

    TEST(model.B[k].rows == model.nstates);
    TEST(model.B[k].cols == model.ninputs);

    TEST(model.f[k].rows == model.nstates);
    TEST(model.f[k].cols == 1);

  }
  TEST(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON-1) < tol);
}

void tiny_InitModelDataMatrix_Ltv_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 2;
  const int NINPUTS = 2;
  const int NHORIZON = 3;
  sfloat A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
  sfloat B_data[] = {1.1, 2.2, 1.3, 2.3};        
  sfloat f_data[] = {2.1, 1.2, 3.3, 1.3};     
  Matrix A[NHORIZON-1];
  Matrix B[NHORIZON-1];     
  Matrix f[NHORIZON-1];        
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);

  sfloat* A_ptr = A_data;
  sfloat* B_ptr = B_data;
  sfloat* f_ptr = f_data;

  for (int k = 0; k < model.nhorizon-1; ++k) {
    A[k] = slap_MatrixFromArray(model.nstates, model.nstates, A_ptr);
    A_ptr += model.nstates * model.nstates;
    B[k] = slap_MatrixFromArray(model.nstates, model.ninputs, B_ptr);
    B_ptr += model.nstates * model.ninputs;
    f[k] = slap_MatrixFromArray(model.nstates, 1, f_ptr);
    f_ptr += model.nstates;      
  }  

  tiny_InitModelDataMatrix_Ltv(&model, A, B, f);

  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(model.A[k].rows == model.nstates);
    TEST(model.A[k].cols == model.nstates);

    TEST(model.B[k].rows == model.nstates);
    TEST(model.B[k].cols == model.ninputs);

    TEST(model.f[k].rows == model.nstates);
    TEST(model.f[k].cols == 1);

  }
  TEST(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON-1) < tol);
}

void tiny_InitModelMemory_Ltv_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 2;
  const int NINPUTS = 2;
  const int NHORIZON = 3;
  sfloat A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
  sfloat B_data[] = {1.1, 2.2, 1.3, 2.3};        
  sfloat f_data[] = {2.1, 1.2, 3.3, 1.3};  
                    
  tiny_LtvModel model;
  tiny_SetModelDims_Ltv(&model, NSTATES, NINPUTS, NHORIZON);
  sfloat data[model.data_size];  // data
  Matrix mats[(NHORIZON-1)*3];  // array of matrices
  tiny_InitModelMemory_Ltv(&model, mats, data);

  // just for testing, use UpdateJacobian istead
  tiny_FillModelMemory_Ltv(&model, A_data, B_data, f_data);  

  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(model.A[k].rows == model.nstates);
    TEST(model.A[k].cols == model.nstates);

    TEST(model.B[k].rows == model.nstates);
    TEST(model.B[k].cols == model.ninputs);

    TEST(model.f[k].rows == model.nstates);
    TEST(model.f[k].cols == 1);

  }
  TEST(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON-1) < tol);
}

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
  tiny_InitModelDataArray_Ltv_Test();
  tiny_InitModelDataMatrix_Ltv_Test();
  tiny_InitModelMemory_Ltv_Test();
  tiny_SetModelJacFunc_Ltv_Test();
  tiny_SetModelNonlinear_Ltv_Test();
  PrintTestResult();
  return TestResult();
}
