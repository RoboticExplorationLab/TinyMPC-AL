#include "tiny_struct.h"
#include "tiny_cost.h"
#include "tiny_dynamics.h"

// A and B are computed before LQR
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtiModel model,
                                         Matrix* Q_temp);

enum slap_ErrorCode tiny_BackwardPassLtv(
    tiny_ProblemData* prob, const tiny_Solver solver,
    const tiny_LtvModel model, Matrix* Q_temp);

// enum slap_ErrorCode tiny_BackwardPassLtvf(
//     tiny_ProblemData* prob, const tiny_Solver solver,
//     const tiny_LtvModel model, Matrix Q_temp);

// enum slap_ErrorCode tiny_ForwardPassLtvf(
//     Matrix* X, Matrix* U, const tiny_ProblemData prob, 
//     const tiny_LtvModel model);