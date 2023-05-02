#include "model_ltv.h"

enum tiny_ErrorCode tiny_SetLtvModelDims(tiny_LtvModel* model, const int nstates,
                                        const int ninputs, const int nhorizon,
                                        const int affine, const sfloat dt) {
  SLAP_ASSERT(nstates > 0 && ninputs > 0 && nstates >= ninputs, SLAP_INVALID_DIMENSION, TINY_SLAP_ERROR, 
  "SetModelDims: nstates (%d) >= ninputs (%d) > 0", nstates, ninputs);  
  SLAP_ASSERT(nhorizon > 0, SLAP_INVALID_DIMENSION, TINY_SLAP_ERROR, 
  "SetModelDims: nhorizon (%d) > 0", nhorizon);   
  model->nstates = nstates;
  model->ninputs = ninputs;
  model->nhorizon = nhorizon;
  model->affine = affine;
  model->dt = dt;
  model->data_size = nstates*(nstates + ninputs)*(nhorizon - 1);
  if (affine) {
    model->data_size += nstates*(nhorizon - 1);
  }
  model->A = TINY_NULL;
  model->B = TINY_NULL;
  model->f = TINY_NULL;
  model->get_jacobians = TINY_NULL;
  model->get_nonl_model = TINY_NULL;  
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitLtvModelDataMatrix(tiny_LtvModel* model, 
    Matrix* A, Matrix* B, Matrix* f) {
  SLAP_ASSERT(A != TINY_NULL && B != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: A and B must not be TINY_NULL");
  model->A = A;
  model->B = B;
  if (model->affine) {
    SLAP_ASSERT(f != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
    "InitModelData: f must not be TINY_NULL");
    if (f) {  
      model->f = f;
    }
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitModelLtvDataArray(tiny_LtvModel* model, Matrix* A, 
    Matrix* B, Matrix* f, sfloat* A_array, sfloat* B_array, sfloat* f_array) {
  SLAP_ASSERT(A != TINY_NULL && B != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: A and B must not be TINY_NULL");
  SLAP_ASSERT(A_array != TINY_NULL && A_array != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: A_array and B_array must not be TINY_NULL");  
  model->A = A;
  model->B = B;
  
  sfloat* A_ptr = A_array;
  sfloat* B_ptr = B_array;

  for (int k = 0; k < model->nhorizon-1; ++k) {
    model->A[k] = slap_MatrixFromArray(model->nstates, model->nstates, A_ptr);
    A_ptr += model->nstates * model->nstates;
    model->B[k] = slap_MatrixFromArray(model->nstates, model->ninputs, B_ptr);
    B_ptr += model->nstates * model->ninputs; 
  }

  if (model->affine) {
  SLAP_ASSERT(f_array != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: f_array must not be TINY_NULL");  
  SLAP_ASSERT(f != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: f must not be TINY_NULL");
    if (f && f_array) {  
      model->f = f;
      sfloat* f_ptr = f_array;
      for (int k = 0; k < model->nhorizon-1; ++k) {
        model->f[k] = slap_MatrixFromArray(model->nstates, 1, f_ptr);
        f_ptr += model->nstates;     
      }
    }
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitLtvModelMemory(tiny_LtvModel* model, Matrix* mats,
    sfloat* data) {
  SLAP_ASSERT(data != TINY_NULL && mats != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelMemory: data and mats must not be TINY_NULL");
  Matrix* mat_ptr = mats;
  model->A = mat_ptr;
  mat_ptr += model->nhorizon - 1;
  model->B = mat_ptr;

  if (model->affine) {
    mat_ptr += model->nhorizon - 1;
    model->f = mat_ptr;
  }
  
  sfloat* ptr = data;
  for (int k = 0; k < model->nhorizon-1; ++k) {
    model->A[k] = slap_MatrixFromArray(model->nstates, model->nstates, ptr);
    ptr += model->nstates * model->nstates;
    model->B[k] = slap_MatrixFromArray(model->nstates, model->ninputs, ptr);
    ptr += model->nstates * model->ninputs;
    if (model->affine) {
      model->f[k] = slap_MatrixFromArray(model->nstates, 1, ptr);
      ptr += model->nstates;      
    }
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_FillLtvModelMemory(tiny_LtvModel* model, sfloat* A_data, 
sfloat* B_data, sfloat* f_data) {
  sfloat* A_ptr = A_data;
  sfloat* B_ptr = B_data;

  for (int k = 0; k < model->nhorizon-1; ++k) {
    model->A[k] = slap_MatrixFromArray(model->nstates, model->nstates, A_ptr);
    A_ptr += model->nstates * model->nstates;
    model->B[k] = slap_MatrixFromArray(model->nstates, model->ninputs, B_ptr);
    B_ptr += model->nstates * model->ninputs;
   
  }

  if (model->affine) {
    sfloat* f_ptr = f_data;
    for (int k = 0; k < model->nhorizon-1; ++k) {
      model->f[k] = slap_MatrixFromArray(model->nstates, 1, f_ptr);
      f_ptr += model->nstates;   
    }
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetLtvModelJacFunc(
    tiny_LtvModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_jacobians != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "SetModelJacFunc: pointer to function must not be TINY_NULL");
  model->get_jacobians = get_jacobians;
  return TINY_NO_ERROR;        
}

enum tiny_ErrorCode tiny_SetLtvModelNonlFunc(
    tiny_LtvModel* model, 
    void (*get_nonl_model)(Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_nonl_model != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "SetModelNonlinear: pointer to function must not be TINY_NULL");
  model->get_nonl_model = get_nonl_model;
  return TINY_NO_ERROR;    
}                      