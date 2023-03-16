#pragma once
#ifdef __cplusplus
extern "C" {
#endif

// NOTE: this is an odd fix to get gcov to run correctly on GitHub Actions:
// https://www.osadl.org/fileadmin/dam/interface/docbook/howtos/coverage.pdf
// void __gcov_flush(void);

#include "tiny_struct.h"
#include "tiny_utils.h"
#include "tiny_cost.h"
#include "tiny_constraint.h"
#include "tiny_dynamics_lti.h"
#include "tiny_mpc_lti.h"
#include "tiny_lqr_lti.h"
#include "tiny_dynamics_ltv.h"
#include "tiny_mpc_ltv.h"
#include "tiny_lqr_ltv.h"

// #if (MODEL == LTI_MODEL)
//   #include "tiny_dynamics_lti.h"
//   #if (CONSTRAINT == CONIC_CONSTRAINT)
//   // #include "tiny_conic_mpc_lti.h"
//     // #include "tiny_conic_constraint.h"
//   #endif
//   #if (CONSTRAINT == LINEAR_CONSTRAINT)
//     #include "tiny_mpc_lti.h"
//     #include "tiny_constraint.h"
//   #endif
//   #if (CONSTRAINT == NO_CONSTRAINT)
//     #include "tiny_lqr_lti.h"
//   #endif
// #endif


// #if (MODEL == LTV_MODEL)
//   #include "tiny_dynamics_ltv.h"
//   #if (CONSTRAINT == CONIC_CONSTRAINT)
//   // #include "tiny_conic_mpc_ltv.h"
//     // #include "tiny_conic_constraint.h"
//   #endif
//   #if (CONSTRAINT == LINEAR_CONSTRAINT)
//     #include "tiny_mpc_ltv.h"
//     #include "tiny_constraint.h"
//   #endif
//   #if (CONSTRAINT == NO_CONSTRAINT)
//     #include "tiny_lqr_ltv.h"
//   #endif
// #endif

#ifdef __cplusplus
}
#endif