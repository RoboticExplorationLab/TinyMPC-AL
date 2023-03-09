#include "tiny_struct.h"

void tiny_InitLinearDiscreteModel(tiny_LinearDiscreteModel* model) {
  *model = (tiny_LinearDiscreteModel){
      .nstates = 0,
      .ninputs = 0,
      .dt = 0.0,
      .A = kNullMat,
      .B = kNullMat,
      .f = kNullMat,
      .x0 = kNullMat,
  };
}

void tiny_InitKnotPoint(tiny_KnotPoint* z) {
  *z = (tiny_KnotPoint){
      .x = kNullMat,
      .u = kNullMat,
      .t = 0.0,
      .dt = 0.0,
  };
}

void tiny_InitSolver(tiny_Solver* solver) {
  *solver = (tiny_Solver){
      .regu = 1e-8,
      .regu_min = 1e-8,
      .regu_max = 1e2,
      .penalty = 1,
      .penalty_max = 1e8,
      .penalty_mul = 10,
      .max_primal_iters = 50,
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
      .Q = kNullMat,
      .R = kNullMat,
      .q = kNullMat,
      .r = kNullMat,
      .Qf = kNullMat,
      .qf = kNullMat,
      .u_max = kNullMat,
      .u_min = kNullMat,
      .x_max = kNullMat,
      .x_min = kNullMat,
      .X_ref = NULL,
      .U_ref = NULL,
      .dt = 0.0,
      .x0 = kNullMat,
      .K = NULL,
      .d = NULL,
      .P = NULL,
      .p = NULL,
      .input_duals = NULL,
      .state_duals = NULL,
      .goal_dual = kNullMat,
  };
}
