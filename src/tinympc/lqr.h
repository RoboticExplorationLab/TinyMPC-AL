#ifndef LQR_H
# define LQR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "cost_lqr.h"
#include "model.h"

enum tiny_ErrorCode tiny_RollOutModelCost(tiny_Workspace* work);
enum tiny_ErrorCode tiny_ForwardPass(tiny_Workspace* work);
enum tiny_ErrorCode tiny_BackwardPass(tiny_Workspace* work);
enum tiny_ErrorCode tiny_SolveLqr(tiny_Workspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef LQR_H