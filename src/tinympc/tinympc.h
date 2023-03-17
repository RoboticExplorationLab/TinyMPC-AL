#pragma once
#ifdef __cplusplus
extern "C" {
#endif

#include "data_struct.h"
#include "utils.h"
#include "cost_lqr.h"
#include "constraint_linear.h"
#include "dynamics_lti.h"
#include "mpc_lti.h"
#include "lqr_lti.h"
#include "dynamics_ltv.h"
#include "mpc_ltv.h"
#include "lqr_ltv.h"

// #if (MODEL == LTI_MODEL)
//   #include "dynamics_lti.h"
//   #if (CONSTRAINT == CONIC_CONSTRAINT)
//   // #include "conic_mpc_lti.h"
//     // #include "conic_constraint.h"
//   #endif
//   #if (CONSTRAINT == LINEAR_CONSTRAINT)
//     #include "mpc_lti.h"
//     #include "constraint.h"
//   #endif
//   #if (CONSTRAINT == NO_CONSTRAINT)
//     #include "lqr_lti.h"
//   #endif
// #endif


// #if (MODEL == LTV_MODEL)
//   #include "dynamics_ltv.h"
//   #if (CONSTRAINT == CONIC_CONSTRAINT)
//   // #include "conic_mpc_ltv.h"
//     // #include "conic_constraint.h"
//   #endif
//   #if (CONSTRAINT == LINEAR_CONSTRAINT)
//     #include "mpc_ltv.h"
//     #include "constraint.h"
//   #endif
//   #if (CONSTRAINT == NO_CONSTRAINT)
//     #include "lqr_ltv.h"
//   #endif
// #endif

#ifdef __cplusplus
}
#endif