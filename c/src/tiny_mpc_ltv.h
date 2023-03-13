#include "tiny_constraint.h"
#include "tiny_cost.h"
#include "tiny_dynamics_ltv.h"
#include "tiny_struct.h"

enum slap_ErrorCode tiny_ConstrainedBackwardPassLtv(
    tiny_ProblemData* prob, const tiny_Solver solver, const tiny_LtvModel model,
    const Matrix* X, const Matrix* U, Matrix* Q_temp, Matrix* ineq_temp);

enum slap_ErrorCode tiny_MpcLtv(Matrix* X, Matrix* U, tiny_ProblemData* prob,
                                tiny_Solver* solver, const tiny_LtvModel model,
                                const int verbose);