#ifndef COST_LQR_H
# define COST_LQR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "utils.h"

enum tiny_ErrorCode tiny_AddStageCost(tiny_Workspace* work, const int k);

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_Workspace* work);

enum tiny_ErrorCode tiny_ExpandStageCost(tiny_Workspace* work, const int k);

enum tiny_ErrorCode tiny_ExpandTerminalCost(tiny_Workspace* work);

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_Workspace* work);

// enum tiny_ErrorCode tiny_SetQ(tiny_Workspace* work, sfloat* Q);

// enum tiny_ErrorCode tiny_SetQf(tiny_Workspace* work, sfloat* Qf);

// enum tiny_ErrorCode tiny_SetR(tiny_Workspace* work, sfloat* R);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef COST_LQR_H