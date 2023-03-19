#include "model_lti.h"

// Calculate data_size as well
enum tiny_ErrorCode tiny_SetModelDims_Lti(tiny_LtiModel* model, 
    const int nstates, const int ninputs) {
  SLAP_ASSERT(nstates > 0 && ninputs > 0 && nstates >= ninputs, SLAP_INVALID_DIMENSION, TINY_SLAP_ERROR, 
  "SetModelDims: nstates (%d) >= ninputs (%d) > 0", nstates, ninputs);
  model->nstates = nstates;
  model->ninputs = ninputs;
  model->data_size = nstates*(nstates + ninputs + 1);
  model->A = kNullMat;
  model->B = kNullMat;
  model->f = kNullMat;
  model->get_jacobians = NULL;
  model->get_nonlinear_dynamics = NULL;
  return TINY_NO_ERROR;
}

// User already has data
enum tiny_ErrorCode tiny_InitModelData_Lti(tiny_LtiModel* model, 
    sfloat* A, sfloat* B, sfloat* f) {
  SLAP_ASSERT(A != NULL && B != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelData: A and B must not be NULL");
  model->A = slap_MatrixFromArray(model->nstates, model->nstates, A);
  model->B = slap_MatrixFromArray(model->nstates, model->ninputs, B);
  model->f = slap_MatrixFromArray(model->nstates, 1, f);
  return TINY_NO_ERROR;
}

// User gives memory and then assign data
enum tiny_ErrorCode tiny_InitModelMemory_Lti(tiny_LtiModel* model, sfloat* data) {
  SLAP_ASSERT(data != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "InitModelMemory: data must not be NULL");
  sfloat* ptr = data;
  model->A = slap_MatrixFromArray(model->nstates, model->nstates, ptr);
  ptr += model->nstates * model->nstates;
  model->B = slap_MatrixFromArray(model->nstates, model->ninputs, ptr);
  ptr += model->nstates * model->ninputs;
  model->f = slap_MatrixFromArray(model->nstates, 1, ptr);  
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetModelJacFunc_Lti(
    tiny_LtiModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_jacobians != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "SetModelJacFunc: pointer to function must not be NULL");
  model->get_jacobians = get_jacobians;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetModelNonlinear_Lti(
    tiny_LtiModel* model, 
    void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_nonlinear_dynamics != NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "SetModelNonlinear: pointer to function must not be NULL");
  model->get_nonlinear_dynamics =*get_nonlinear_dynamics;
  return TINY_NO_ERROR;
}