#pragma once

#include "utils.h"

typedef struct {
  int nstates;
  int ninputs;
  int nhorizon;
  Matrix* A;
  Matrix* B;
  Matrix* f;
  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix);
  int data_size;
} tiny_LtvModel;

enum tiny_ErrorCode tiny_SetModelDims_Ltv(tiny_LtvModel* model, const int nstates,
                                        const int ninputs, const int nhorizon);

enum tiny_ErrorCode tiny_InitModelData_Ltv(tiny_LtvModel* model, 
    sfloat* A, sfloat* B, sfloat* f);

enum tiny_ErrorCode tiny_InitModelMemory_Ltv(tiny_LtvModel* model, sfloat* data);

enum tiny_ErrorCode tiny_SetModelJacFunc_Ltv(
    tiny_LtvModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix));

enum tiny_ErrorCode tiny_SetModelNonlinear_Ltv(
    tiny_LtvModel* model, 
    void (*get_nonlinear_dynamics)(Matrix*, Matrix*, const Matrix, 
    const Matrix));                      