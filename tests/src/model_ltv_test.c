#include "simpletest.h"
#include "test_utils.h"
#include "tinympc/model_ltv.h"
#include "tinympc/utils.h"
#include "bicycle_5d.h"

void tiny_SetLtvModelDims_Test() {
  // Test good 
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_LtvModel model;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1)*(NHORIZON-1);
  TEST(tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1) == TINY_NO_ERROR &&
       true_size == model.data_size);

  tiny_LtvModel model_not_affine;
  TEST(tiny_SetLtvModelDims(&model_not_affine, NSTATES, NINPUTS, NHORIZON, 0, 0.1) == TINY_NO_ERROR);
  TEST(model_not_affine.dt == (sfloat)(0.1) && model_not_affine.affine == 0 && 
       true_size - NSTATES*(NHORIZON-1) == model_not_affine.data_size);
  // Test bad (there are more bad cases)
  // NSTATES = 1;
  // NINPUTS = 0;
  // int NHORIZON = 0;
  // TEST(tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1) == TINY_SLAP_ERROR);
}

void tiny_InitModelLtvDataArray_Test() {
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
  tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);
  tiny_InitModelLtvDataArray(&model, A, B, f, A_data, B_data, f_data);
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

  tiny_LtvModel model_not_affine;
  tiny_SetLtvModelDims(&model_not_affine, NSTATES, NINPUTS, NHORIZON, 0, 0.1);
  tiny_InitModelLtvDataArray(&model_not_affine, A, B, TINY_NULL, A_data, B_data, TINY_NULL);
  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(model_not_affine.A[k].rows == model_not_affine.nstates);
    TEST(model_not_affine.A[k].cols == model_not_affine.nstates);

    TEST(model_not_affine.B[k].rows == model_not_affine.nstates);
    TEST(model_not_affine.B[k].cols == model_not_affine.ninputs);
  }
  TEST(SumOfSquaredErrorMatrices(A_data, model_not_affine.A, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(B_data, model_not_affine.B, NHORIZON-1) < tol);
  TEST(model_not_affine.f == TINY_NULL);
}

void tiny_InitLtvModelDataMatrix_Test() {
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
  tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);

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

  tiny_InitLtvModelDataMatrix(&model, A, B, f);

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

  tiny_LtvModel model_not_affine;
  tiny_SetLtvModelDims(&model_not_affine, NSTATES, NINPUTS, NHORIZON, 0, 0.1);
  tiny_InitLtvModelDataMatrix(&model_not_affine, A, B, f);
  TEST(SumOfSquaredErrorMatrices(A_data, model_not_affine.A, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(B_data, model_not_affine.B, NHORIZON-1) < tol);
  TEST(model_not_affine.f == TINY_NULL);
}

void tiny_InitLtvModelMemory_Test() {
  const sfloat tol = 1e-8;
  const int NSTATES = 2;
  const int NINPUTS = 2;
  const int NHORIZON = 3;
  sfloat A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
  sfloat B_data[] = {1.1, 2.2, 1.3, 2.3};        
  sfloat f_data[] = {2.1, 1.2, 3.3, 1.3};  
                    
  tiny_LtvModel model;
  tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);
  sfloat data[model.data_size];  // data
  Matrix mats[(NHORIZON-1)*3];  // array of matrices A, B, f
  tiny_InitLtvModelMemory(&model, mats, data);

  // just for testing, use UpdateJacobian istead
  tiny_FillLtvModelMemory(&model, A_data, B_data, f_data);  

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

  tiny_LtvModel model_not_affine;
  tiny_SetLtvModelDims(&model_not_affine, NSTATES, NINPUTS, NHORIZON, 0, 0.1);
  sfloat data2[model_not_affine.data_size];  // data
  Matrix mats2[(NHORIZON-1)*2];  // array of matrices A and B
  tiny_InitLtvModelMemory(&model_not_affine, mats2, data2);
  tiny_FillLtvModelMemory(&model_not_affine, A_data, B_data, f_data);  
  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(model_not_affine.A[k].rows == model.nstates);
    TEST(model_not_affine.A[k].cols == model.nstates);

    TEST(model_not_affine.B[k].rows == model.nstates);
    TEST(model_not_affine.B[k].cols == model.ninputs);
  }
  TEST(SumOfSquaredErrorMatrices(A_data, model_not_affine.A, NHORIZON-1) < tol);
  TEST(SumOfSquaredErrorMatrices(B_data, model_not_affine.B, NHORIZON-1) < tol);
  TEST(model_not_affine.f == TINY_NULL);
}

void tiny_SetLtvModelJacFunc_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_LtvModel model;
  tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);
  tiny_SetLtvModelJacFunc(&model, tiny_Bicycle5dGetJacobians);
  TEST(model.get_jacobians == tiny_Bicycle5dGetJacobians);
}

void tiny_SetLtvModelNonlFunc_Test() {
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_LtvModel model;
  tiny_SetLtvModelDims(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);
  tiny_SetLtvModelNonlFunc(&model, tiny_Bicycle5dNonlinearDynamics);
  TEST(model.get_nonl_model == tiny_Bicycle5dNonlinearDynamics);
}

int main() {
  tiny_SetLtvModelDims_Test();
  tiny_InitModelLtvDataArray_Test();
  tiny_InitLtvModelDataMatrix_Test();
  tiny_InitLtvModelMemory_Test();
  tiny_SetLtvModelJacFunc_Test();
  tiny_SetLtvModelNonlFunc_Test();
  PrintTestResult();
  return TestResult();
}