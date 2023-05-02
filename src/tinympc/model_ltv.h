#ifndef MODEL_LTV_H
# define MODEL_LTV_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

enum tiny_ErrorCode tiny_SetLtvModelDims(tiny_LtvModel* model, const int nstates,
                                        const int ninputs, const int nhorizon,
                                        const int affine, const sfloat dt);

// User provides array of slap matrices
enum tiny_ErrorCode tiny_InitLtvModelDataMatrix(tiny_LtvModel* model, 
    Matrix* A, Matrix* B, Matrix* f);

// User provides matrix as column-major array
enum tiny_ErrorCode tiny_InitModelLtvDataArray(tiny_LtvModel* model, Matrix* A, 
    Matrix* B, Matrix* f, sfloat* A_array, sfloat* B_array, sfloat* f_array);

enum tiny_ErrorCode tiny_InitLtvModelMemory(tiny_LtvModel* model, Matrix* mats,
    sfloat* data);

// Used after tiny_InitLtvModelMemory and before tiny_UpdateLtvModelJac
enum tiny_ErrorCode tiny_FillLtvModelMemory(tiny_LtvModel* model, sfloat* A_data, 
sfloat* B_data, sfloat* f_data);

enum tiny_ErrorCode tiny_SetLtvModelJacFunc(
    tiny_LtvModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix));

enum tiny_ErrorCode tiny_SetLtvModelNonlFunc(
    tiny_LtvModel* model, 
    void (*get_nonl_model)(Matrix*, const Matrix, const Matrix));


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef MODEL_LTV_H