#pragma once

#include "utils.h"

typedef struct {
  int nstates;
  int ninputs;
  Matrix A;
  Matrix B;
  Matrix f;
  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix);
  int data_size;  ///< number of doubles need to store the data
} tiny_LtiModel;

enum tiny_ErrorCode tiny_SetModelDims_Lti(tiny_LtiModel* model, 
    const int nstates, const int ninputs);

enum tiny_ErrorCode tiny_InitModelData_Lti(tiny_LtiModel* model, 
    sfloat* A, sfloat* B, sfloat* f);

enum tiny_ErrorCode tiny_InitModelMemory_Lti(tiny_LtiModel* model, sfloat* data);

enum tiny_ErrorCode tiny_SetModelJacFunc_Lti(
    tiny_LtiModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix));

enum tiny_ErrorCode tiny_SetModelNonlinear_Lti(
    tiny_LtiModel* model, 
    void (*get_nonlinear_dynamics)(Matrix*, Matrix*, const Matrix, 
    const Matrix));
