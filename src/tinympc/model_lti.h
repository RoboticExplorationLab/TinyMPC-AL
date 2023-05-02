#ifndef MODEL_LTI_H
# define MODEL_LTI_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"

enum tiny_ErrorCode tiny_SetLtiModelDims(tiny_LtiModel* model, 
                                         const int nstates, 
                                         const int ninputs,
                                         const int affine,
                                         const sfloat dt);

enum tiny_ErrorCode tiny_InitLtiModelData(tiny_LtiModel* model, 
                                          sfloat* A, 
                                          sfloat* B, 
                                          sfloat* f);

enum tiny_ErrorCode tiny_InitLtiModelMemory(tiny_LtiModel* model, 
                                            sfloat* data);

enum tiny_ErrorCode tiny_SetLtiModelJacFunc(
    tiny_LtiModel* model, 
    void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix));

enum tiny_ErrorCode tiny_SetLtiModelNonlFunc(
    tiny_LtiModel* model, 
    void (*get_nonl_model)(Matrix*, const Matrix, const Matrix));

// enum tiny_ErrorCode tiny_EvalLtiModel(Matrix* xn, const Matrix x, const Matrix u,
//                       const tiny_LtiModel model);

// enum slap_ErrorCode tiny_RollOutLtiModel(tiny_Workspace* work,
//                                          const tiny_LtiModel model);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef MODEL_LTI_H
