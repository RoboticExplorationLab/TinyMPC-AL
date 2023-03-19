#include "model_ltv.h"

enum tiny_ErrorCode tiny_SetModelDims_Ltv(tiny_LtvModel* model, const int nstates,
                                        const int ninputs, const int nhorizon) {
  SLAP_ASSERT(nstates > 0 && ninputs > 0 && nstates >= ninputs, SLAP_INVALID_DIMENSION, SLAP_INVALID_DIMENSION, 
  "SetModelDims: nstates (%d) >= ninputs (%d) > 0", nstates, ninputs);  
  SLAP_ASSERT(nhorizon > 0, SLAP_INVALID_DIMENSION, SLAP_INVALID_DIMENSION, 
  "SetModelDims: nhorizon (%d) > 0", nhorizon);   
  model->nstates = nstates;
  model->ninputs = ninputs;
  model->nhorizon = nhorizon;
  model->data_size = nstates*(nstates + ninputs + 1)*(nhorizon - 1);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitModelData_Ltv(tiny_LtvModel* model, 
    sfloat* A, sfloat* B, sfloat* f) {
  SLAP_ASSERT(A != NULL && B != NULL, SLAP_BAD_POINTER, SLAP_BAD_POINTER,
  "InitModelData: A or B must not be NULL");
  sfloat* Aptr = A;
  sfloat* Bptr = B;
  sfloat* fptr = f;
  for (int k = 0; k < model->nhorizon - 1; ++k) {
    model->A[k] = slap_MatrixFromArray(model->nstates, model->nstates, Aptr);
    Aptr += model->nstates * model->nstates;
    model->B[k] = slap_MatrixFromArray(model->nstates, model->ninputs, Bptr);
    Bptr += model->nstates * model->ninputs;
    model->f[k] = slap_MatrixFromArray(model->nstates, 1, fptr);
    fptr += model->nstates;
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitModelMemory_Ltv(tiny_LtvModel* model, sfloat* data) {
  SLAP_ASSERT(data != NULL, SLAP_BAD_POINTER, SLAP_BAD_POINTER,
  "InitModelMemory: data must not be NULL");
  sfloat* ptr = data;
  for (int k = 0; k < model->nhorizon - 1; ++k) {
    model->A[k] = slap_MatrixFromArray(model->nstates, model->nstates, ptr);
    ptr += model->nstates * model->nstates;
    model->B[k] = slap_MatrixFromArray(model->nstates, model->ninputs, ptr);
    ptr += model->nstates * model->ninputs;
    model->f[k] = slap_MatrixFromArray(model->nstates, 1, ptr);
    ptr += model->nstates;      
    return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetModelJacFunc_Ltv(
    tiny_LtvModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_jacobians != NULL, SLAP_BAD_POINTER, SLAP_BAD_POINTER,
  "SetModelJacFunc: pointer to function must not be NULL");
  model->get_jacobians = get_jacobians;
  return TINY_NO_ERROR;        
}

enum tiny_ErrorCode tiny_SetModelNonlinear_Ltv(
    tiny_LtvModel* model, 
    void (*get_nonlinear_dynamics)(Matrix*, Matrix*, const Matrix, 
    const Matrix)) {
  SLAP_ASSERT(get_nonlinear_dynamics != NULL, SLAP_BAD_POINTER, SLAP_BAD_POINTER,
  "SetModelNonlinear: pointer to function must not be NULL");
  model->get_nonlinear_dynamics = get_nonlinear_dynamics;
  return TINY_NO_ERROR;    
}                      