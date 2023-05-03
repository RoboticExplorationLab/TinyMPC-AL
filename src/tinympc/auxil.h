#ifndef AUXIL_H
# define AUXIL_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"

enum tiny_ErrorCode tiny_InitSettings(tiny_Settings* stgs);

enum tiny_ErrorCode tiny_InitData(tiny_Workspace* work);

enum tiny_ErrorCode tiny_InitSolution(tiny_Workspace* work);

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
