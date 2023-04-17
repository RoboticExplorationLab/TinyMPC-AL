#include "cost_lqr.h"
#include "data_struct.h"
#include "dynamics_ltv.h"

enum slap_ErrorCode tiny_BackwardPassLtv(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtvModel model,
                                         Matrix* Q_temp);

// TODO: Allow online A, B computation? (require less memory but more time)
// enum slap_ErrorCode tiny_BackwardPassLtvf(
//     tiny_ProblemData* prob, const tiny_Solver solver,
//     const tiny_LtvModel model, Matrix* Q_temp);