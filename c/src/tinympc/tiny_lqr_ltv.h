#include "tiny_cost.h"
#include "tiny_dynamics_ltv.h"
#include "tiny_struct.h"

enum slap_ErrorCode tiny_BackwardPassLtv(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtvModel model,
                                         Matrix* Q_temp);

// TODO: Allow online A, B computation? (require less memory but more time)
// enum slap_ErrorCode tiny_BackwardPassLtvf(
//     tiny_ProblemData* prob, const tiny_Solver solver,
//     const tiny_LtvModel model, Matrix* Q_temp);