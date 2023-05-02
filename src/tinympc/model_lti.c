#include "model_lti.h"

// Calculate data_size as well
enum tiny_ErrorCode tiny_SetLtiModelDims(tiny_LtiModel* model, 
    const int nstates, const int ninputs, const int affine, const sfloat dt) 
{
  SLAP_ASSERT(nstates > 0 && ninputs > 0 && nstates >= ninputs, SLAP_INVALID_DIMENSION, TINY_SLAP_ERROR, 
  "tiny_SetLtiModelDims: nstates (%d) >= ninputs (%d) > 0", nstates, ninputs);
  model->nstates = nstates;
  model->ninputs = ninputs;
  model->affine  = affine;
  model->dt      = dt;
  model->data_size = nstates*(nstates + ninputs);
  if (model->affine) {
    model->data_size += nstates;
  }
  model->A = TINY_NULL_MAT;
  model->B = TINY_NULL_MAT;
  model->f = TINY_NULL_MAT;
  model->get_jacobians = TINY_NULL;
  model->get_nonl_model = TINY_NULL;
  return TINY_NO_ERROR;
}

// User already has data
enum tiny_ErrorCode tiny_InitLtiModelData(tiny_LtiModel* model, 
    sfloat* A, sfloat* B, sfloat* f) {
  SLAP_ASSERT(A != TINY_NULL && B != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitLtiModelData: A and B must not be TINY_NULL");
  model->A = slap_MatrixFromArray(model->nstates, model->nstates, A);
  model->B = slap_MatrixFromArray(model->nstates, model->ninputs, B);
  if (model->affine) {
    SLAP_ASSERT(f != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
    "tiny_InitLtiModelData: f must not be TINY_NULL");
    if (f) {
        model->f = slap_MatrixFromArray(model->nstates, 1, f);
    }
  }
  return TINY_NO_ERROR;
}

// User gives memory and then assign data
enum tiny_ErrorCode tiny_InitLtiModelMemory(tiny_LtiModel* model, sfloat* data) {
  SLAP_ASSERT(data != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_InitLtiModelMemory_Lti: data must not be TINY_NULL");
  sfloat* ptr = data;
  model->A = slap_MatrixFromArray(model->nstates, model->nstates, ptr);
  ptr += model->nstates * model->nstates;
  model->B = slap_MatrixFromArray(model->nstates, model->ninputs, ptr);
  ptr += model->nstates * model->ninputs;
  if (model->affine) {
    model->f = slap_MatrixFromArray(model->nstates, 1, ptr);  
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetLtiModelJacFunc(
    tiny_LtiModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_jacobians != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_SetLtiModelJacFunc: pointer to function must not be TINY_NULL");
  model->get_jacobians = get_jacobians;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetLtiModelNonlFunc(
    tiny_LtiModel* model, 
    void (*get_nonl_model)(Matrix*, const Matrix, const Matrix)) {
  SLAP_ASSERT(get_nonl_model != TINY_NULL, SLAP_BAD_POINTER, TINY_SLAP_ERROR,
  "tiny_SetLtiModelNonlinear: pointer to function must not be TINY_NULL");
  model->get_nonl_model = get_nonl_model;
  return TINY_NO_ERROR;
}

// enum slap_ErrorCode tiny_EvalLtiModel(Matrix* xn, const Matrix x, const Matrix u,
//                       const tiny_LtiModel model) {
//   slap_Copy(*xn, model.f);
//   slap_MatMulAdd(*xn, model.A, x, 1, 1);  // x[k+1] += A * x[k]
//   slap_MatMulAdd(*xn, model.B, u, 1, 1);  // x[k+1] += B * u[k]
//   return SLAP_NO_ERROR;
// }

// enum slap_ErrorCode tiny_RollOutLtiModel(tiny_Workspace* work,
//                                          const tiny_LtiModel model) {
//   int N = work->data->nhorizon;
//   for (int k = 0; k < N - 1; ++k) {
//     // delta_x and delta_u over previous X, U
//     // Control input: u = uf - d - K*(x - xf)
//     slap_Copy(work->soln->U[k], work->soln->d[k]);                     // u[k] = d[k]
//     slap_MatMulAdd(work->soln->U[k], work->soln->.K[k], work->soln->X[k], -1, -1);  // u[k] += K[k] * x[k]
//     // Next state: x = A*x + B*u + f
//     tiny_EvalLtiModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], model);
//   }
//   return SLAP_NO_ERROR;
// }