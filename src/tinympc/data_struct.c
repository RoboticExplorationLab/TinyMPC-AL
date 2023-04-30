#include "data_struct.h"

void tiny_InitLtiModel(tiny_LtiModel* model) {
  *model = (tiny_LtiModel){
      .nstates = 0,
      .ninputs = 0,
      .dt = 0.0,
      .A = TINY_NULL_MAT,
      .B = TINY_NULL_MAT,
      .f = TINY_NULL_MAT,
      .x0 = TINY_NULL_MAT,
  };
}

void tiny_InitLtvModel(tiny_LtvModel* model) {
  *model = (tiny_LtvModel){
      .nstates = 0,
      .ninputs = 0,
      .dt = 0.0,
      .A = NULL,
      .B = NULL,
      .f = NULL,
      .x0 = TINY_NULL_MAT,
      .get_jacobians = NULL,
      .get_nonl_model = NULL,
  };
}

void tiny_InitKnotPoint(tiny_KnotPoint* z) {
  *z = (tiny_KnotPoint){
      .x = TINY_NULL_MAT,
      .u = TINY_NULL_MAT,
      .t = 0.0,
      .dt = 0.0,
  };
}

void tiny_InitSettings(tiny_Settings* solver) {
  *solver = (tiny_Settings){
      .reg = 1e-8,
      .reg_min = 1e-8,
      .reg_max = 1e2,
      .penalty = 1,
      .penalty_max = 1e8,
      .penalty_mul = 10,
      .max_outer_iters = 50,
      .max_search_iters = 5,
      .riccati_tol = 1e-1,
      .cstr_tol = 1e-4,
  };
}

void tiny_InitProblemData(tiny_ProblemData* prob) {
  *prob = (tiny_ProblemData){
      .nstates = 0,
      .ninputs = 0,
      .nhorizon = 0,
      .ncstr_states = 0,
      .ncstr_inputs = 0,
      .ncstr_goal = 0,
      .Q = TINY_NULL_MAT,
      .R = TINY_NULL_MAT,
      .q = TINY_NULL_MAT,
      .r = TINY_NULL_MAT,
      .Qf = TINY_NULL_MAT,
      .qf = TINY_NULL_MAT,
      .X_ref = NULL,
      .U_ref = NULL,
      .dt = 0.0,
      .x0 = TINY_NULL_MAT,
      .K = NULL,
      .d = NULL,
      .P = NULL,
      .p = NULL,
      .YU = NULL,
      .YX = NULL,
      .YG = TINY_NULL_MAT,
      .Acx = TINY_NULL_MAT,
      .bcx = TINY_NULL_MAT,
      .Acu = TINY_NULL_MAT,
      .bcu = TINY_NULL_MAT,
  };
}
