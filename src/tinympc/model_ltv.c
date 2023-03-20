#include "model_ltv.h"

enum tiny_ErrorCode tiny_SetModelDims_Ltv(tiny_LtvModel* model, const int nstates,
                                        const int ninputs, const int nhorizon) {
  SLAP_ASSERT(nstates > 0 && ninputs > 0 && nstates >= ninputs, SLAP_INVALID_DIMENSION, TINY_SLAP_ERROR, 
  "SetModelDims: nstates (%d) >= ninputs (%d) > 0", nstates, ninputs);  
  SLAP_ASSERT(nhorizon > 0, SLAP_INVALID_DIMENSION, TINY_SLAP_ERROR, 
  "SetModelDims: nhorizon (%d) > 0", nhorizon);   
  model->nstates = nstates;
  model->ninputs = ninputs;
  model->nhorizon = nhorizon;
  model->data_size = nstates*(nstates + ninputs + 1)*(nhorizon - 1);
  model->A = NULL;
  model->B = NULL;
  model->f = NULL;
  model->get_jacobians = NULL;
  model->get_nonlinear_dynamics = NULL;  
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitModelData_Ltv(tiny_LtvModel* model, 
    Matrix* A, Matrix* B, Matrix* f) {
  SLAP_ASSERT(A != NULL && B != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: A and B must not be NULL");
  model->A = A;
  model->B = B;
  model->f = f;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitModelMemory_Ltv(tiny_LtvModel* model, Matrix* mats,
    sfloat* data) {
  SLAP_ASSERT(data != NULL && mats != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelMemory: data and mats must not be NULL");
  Matrix* mat_ptr = mats;
  model->A = mat_ptr;
  mat_ptr += model->nhorizon - 1;
  model->B = mat_ptr;
  mat_ptr += model->nhorizon - 1;
  model->f = mat_ptr;
  
  sfloat* ptr = data;
  for (int k = 0; k < model->nhorizon - 1; ++k) {
    model->A[k] = slap_MatrixFromArray(model->nstates, model->nstates, ptr);
    ptr += model->nstates * model->nstates;
    model->B[k] = slap_MatrixFromArray(model->nstates, model->ninputs, ptr);
    ptr += model->nstates * model->ninputs;
    model->f[k] = slap_MatrixFromArray(model->nstates, 1, ptr);
    ptr += model->nstates;      
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetModelJacFunc_Ltv(
    tiny_LtvModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_jacobians != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "SetModelJacFunc: pointer to function must not be NULL");
  model->get_jacobians = get_jacobians;
  return TINY_NO_ERROR;        
}

enum tiny_ErrorCode tiny_SetModelNonlinear_Ltv(
    tiny_LtvModel* model, 
    void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_nonlinear_dynamics != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "SetModelNonlinear: pointer to function must not be NULL");
  model->get_nonlinear_dynamics = get_nonlinear_dynamics;
  return TINY_NO_ERROR;    
}                      