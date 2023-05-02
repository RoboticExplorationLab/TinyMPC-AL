#ifndef MODEL_H
# define MODEL_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

enum tiny_ErrorCode tiny_InitModel(tiny_Model* model, const int nstates,
                                   const int ninputs, const int nhorizon,
                                   const int ltv, const int affine, 
                                   const sfloat dt);

// User provides array of slap matrices
enum tiny_ErrorCode tiny_InitModelDataMatrix(tiny_Model* model, 
    Matrix* A, Matrix* B, Matrix* f);

// User provides matrix as column-major array
enum tiny_ErrorCode tiny_InitModelDataArray(tiny_Model* model, Matrix* A, 
    Matrix* B, Matrix* f, sfloat* A_array, sfloat* B_array, sfloat* f_array);

enum tiny_ErrorCode tiny_InitModelMemory(tiny_Model* model, Matrix* mats,
    sfloat* data);

// Used after tiny_InitLtvModelMemory and before tiny_UpdateLtvModelJac
enum tiny_ErrorCode tiny_FillModelMemory(tiny_Model* model, sfloat* A_data, 
sfloat* B_data, sfloat* f_data);

enum tiny_ErrorCode tiny_SetModelJacFunc(
    tiny_Model* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix));

enum tiny_ErrorCode tiny_SetModelNonlFunc(
    tiny_Model* model, 
    void (*get_nonl_model)(Matrix*, const Matrix, const Matrix));

// k = 0 to use LTI model
enum tiny_ErrorCode tiny_EvalModel(Matrix* xn, const Matrix x, const Matrix u,
                                   tiny_Model* model, const int k);

enum tiny_ErrorCode tiny_RollOutClosedLoop(tiny_Workspace* work);

enum tiny_ErrorCode tiny_RollOutOpenLoop(tiny_Workspace* work);

enum tiny_ErrorCode tiny_UpdateModelJac(tiny_Workspace* work);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef MODEL_H