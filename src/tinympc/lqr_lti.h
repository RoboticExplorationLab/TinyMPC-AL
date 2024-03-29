#pragma once
#include "cost_lqr.h"
#include "dynamics_lti.h"
#include "data_struct.h"

// A and B are computed before LQR
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtiModel model,
                                         Matrix* Q_temp);
