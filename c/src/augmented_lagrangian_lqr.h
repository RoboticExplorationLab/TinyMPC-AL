#pragma once

#include "slap/slap.h"

typedef struct tiny_ProblemData {
  int nstates;  
  int ninputs;
  int nhorizon;
  int ncstr_states;
  int ncstr_inputs;
  Matrix Q;
  Matrix R;
  Matrix q;
  Matrix r;
  Matrix Qf;
  Matrix u_max;
  Matrix u_min;
  Matrix x_max;
  Matrix x_min;
  Matrix* x_ref;
  Matrix* u_ref;
  double dt;
  Matrix x0;
  Matrix* K;
  Matrix* d;
  Matrix* P;
  Matrix* p;
} tiny_ProblemData;

const Matrix kNullMat = {
  .rows = 0,
  .cols = 0,
  .sy = 0,
  .is_transposed = 0,
  .data = NULL,
  .mattype = slap_DENSE,
};

const tiny_ProblemData kDefaultProblemData = {
  .nstates = 0,
  .ninputs = 0,
  .nhorizon = 0,
  .ncstr_states = 0,
  .ncstr_inputs = 0,
  .Q = kNullMat,
  .R = kNullMat,
  .q = kNullMat,
  .r = kNullMat,
  .Qf = kNullMat,
  .u_max = kNullMat,
  .u_min = kNullMat,
  .x_max = kNullMat,
  .x_min = kNullMat,
  .x_ref = NULL,
  .u_ref = NULL,
  .dt = 0.0,
  .x0 = kNullMat,
  .K = NULL,
  .d = NULL,
  .P = NULL,
  .p = NULL,  
};

typedef struct {
  Matrix x;  ///< state vector
  Matrix u;  ///< control input vector
  double t;  ///< time
  double dt;  ///< time step
} tiny_KnotPoint;

const tiny_KnotPoint kKnotPoint = {
  .x = kNullMat,
  .u = kNullMat,
  .t = 0.0,
  .dt = 0.0,
};

typedef struct {
  double regu;
  Matrix* input_duals;
  Matrix* state_duals;
  double penalty_min;
  double penalty_max;
  double penalty_mul;
} tiny_Solver;

tiny_Solver kSolver = {
  .regu = 0.0,
  .input_duals = NULL,
  .state_duals = NULL,
  .penalty_min = 0.0,
  .penalty_max = 0.0,
  .penalty_mul = 0.0,
};

typedef struct {
  int nstates;  
  int ninputs;
  double dt;
  Matrix A;
  Matrix B;
  Matrix f;
} tiny_LinearDiscreteModel;

tiny_LinearDiscreteModel kLinearDiscreteModel = {
  .nstates = 0,
  .ninputs = 0,
  .dt = 0.0,
  .A = kNullMat,
  .B = kNullMat,
  .f = kNullMat,
};

enum slap_ErrorCode tiny_BackwardPass(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    const tiny_KnotPoint* Z, const tiny_Solver solver, Matrix S_temp);

enum slap_ErrorCode tiny_ForwardPass(
    const tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    tiny_KnotPoint* Z, const tiny_Solver solver);

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model,
    tiny_KnotPoint* Z, const int verbose);

Matrix slap_IneqInputs(
    const tiny_ProblemData prob, const tiny_KnotPoint* Z);

Matrix slap_IneqInputsJacobian(
    const tiny_ProblemData prob, const tiny_KnotPoint* Z);

Matrix slap_ActiveIneqMask(
    const tiny_Solver solver, const Matrix ineq);    