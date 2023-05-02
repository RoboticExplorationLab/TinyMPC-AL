#pragma once
#include "cost_lqr.h"
#include "types.h"
#include "model.h"

// A and B are computed before LQR
enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Settings solver,
                                         const tiny_LtiModel model,
                                         Matrix* Q_temp);
