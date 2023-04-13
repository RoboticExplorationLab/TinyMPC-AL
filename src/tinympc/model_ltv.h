#pragma once

#include "errors.h"
#include "utils.h"

// for a horizon of N x(0)->x(N-1), need N-1 matrices
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

// User provides array of slap matrices
enum tiny_ErrorCode tiny_InitModelDataMatrix_Ltv(tiny_LtvModel* model, 
    Matrix* A, Matrix* B, Matrix* f);

// User provides matrix as column-major array
enum tiny_ErrorCode tiny_InitModelDataArray_Ltv(tiny_LtvModel* model, Matrix* A, 
    Matrix* B, Matrix* f, sfloat* A_array, sfloat* B_array, sfloat* f_array);

enum tiny_ErrorCode tiny_InitModelMemory_Ltv(tiny_LtvModel* model, Matrix* mats,
    sfloat* data);

enum tiny_ErrorCode tiny_FillModelMemory_Ltv(tiny_LtvModel* model, sfloat* A_data, 
sfloat* B_data, sfloat* f_data);

enum tiny_ErrorCode tiny_SetModelJacFunc_Ltv(
    tiny_LtvModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix));

enum tiny_ErrorCode tiny_SetModelNonlinear_Ltv(
    tiny_LtvModel* model, 
    void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix));