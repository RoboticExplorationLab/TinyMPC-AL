#pragma once
#include "tiny_struct.h"
#include "tiny_cost.h"
#include "tiny_dynamics_lti.h"

// A and B are computed before LQR
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LtiModel model,
                                         Matrix* Q_temp);


