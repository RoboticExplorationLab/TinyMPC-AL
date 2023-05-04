#ifndef AUXIL_H
# define AUXIL_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "utils.h"

enum tiny_ErrorCode tiny_InitSettings(tiny_Settings* stgs);

enum tiny_ErrorCode tiny_SetUnconstrained(tiny_Settings* stgs);

enum tiny_ErrorCode tiny_InitData(tiny_Workspace* work);

enum tiny_ErrorCode tiny_InitDataFromMatrix(tiny_Workspace* work, Matrix x0,
                                            Matrix Q, Matrix R, Matrix Qf,
                                            Matrix* q, Matrix* r, Matrix qf,
                                            Matrix* X_ref, Matrix* U_ref,
                                            Matrix Acx, Matrix bcx,
                                            Matrix Acu, Matrix bcu);

// enum tiny_ErrorCode tiny_InitDataFromArray(tiny_Workspace* work, sfloat* x0,
//                                             sfloat* Q, sfloat* R, sfloat* Qf,
//                                             sfloat* q, sfloat** r, sfloat* qf,
//                                             sfloat* X_ref, sfloat* U_ref,
//                                             sfloat* Acx, sfloat* bcx,
//                                             sfloat* Acu, sfloat* bcu);

enum tiny_ErrorCode tiny_InitDataQuadCostFromArray(tiny_Workspace* work, 
sfloat* Q_data, 
sfloat* R_data, 
sfloat* Qf_data);

enum tiny_ErrorCode tiny_InitDataLinearCostFromArray(tiny_Workspace* work, 
Matrix* q, Matrix* r, 
sfloat* q_data, sfloat* r_data, sfloat* qf_data);

enum tiny_ErrorCode tiny_InitSolution(tiny_Workspace* work);

enum tiny_ErrorCode tiny_InitSolutionFromMatrix(tiny_Workspace* work, 
                                               Matrix* X, Matrix* U,
                                               Matrix* P, Matrix* p,
                                               Matrix* K, Matrix* d,
                                               Matrix* YX, Matrix* YU,
                                               Matrix YG);

enum tiny_ErrorCode tiny_InitSolnTrajFromArray(tiny_Workspace* work,
Matrix* X, Matrix* U,
sfloat* X_data, sfloat* U_data);

enum tiny_ErrorCode tiny_InitSolnDualsFromArray(tiny_Workspace* work,
Matrix* YX, Matrix* YU,
sfloat* YX_data, sfloat* YU_data, sfloat* YG_data);

enum tiny_ErrorCode tiny_InitSolnGainsFromArray(tiny_Workspace* work, 
Matrix* K, Matrix* d, Matrix* P, Matrix* p, 
sfloat* K_data, sfloat* d_data, sfloat* P_data, sfloat* p_data);

enum tiny_ErrorCode tiny_InitWorkspace(tiny_Workspace* work,
                                       tiny_Info* info,
                                       tiny_Model* model,
                                       tiny_Data* data,
                                       tiny_Solution* soln,
                                       tiny_Settings* stgs);

enum tiny_ErrorCode tiny_InitTempData(tiny_Workspace* work, sfloat* temp_data);

enum tiny_ErrorCode tiny_ResetInfo(tiny_Workspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef AUXIL_H
