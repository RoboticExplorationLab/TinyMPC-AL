#pragma once

#include "constraint_linear.h"
#include "cost_lqr.h"
#include "data_struct.h"
#include "dynamics_lti.h"
#include "utils.h"

enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_Solver solver, const tiny_LtiModel model,
    const Matrix* X, const Matrix* U, Matrix* Q_temp, Matrix* ineq_temp);

enum slap_ErrorCode tiny_MpcLti(Matrix* X, Matrix* U, tiny_ProblemData* prob,
                                tiny_Solver* solver, const tiny_LtiModel model,
                                const int verbose, sfloat* temp_data);
