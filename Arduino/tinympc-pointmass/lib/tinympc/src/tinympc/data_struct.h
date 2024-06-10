#pragma once

#include "slap/slap.h"

#define kNullMat  \
  ((Matrix){      \
      0,          \
      0,          \
      0,          \
      0,          \
      NULL,       \
      slap_DENSE, \
  })

typedef struct {
  int nstates;
  int ninputs;
  sfloat dt;
  Matrix A;
  Matrix B;
  Matrix f;
  Matrix x0;
  // int data_size;  ///< number of sfloats need to store the data
} tiny_LtiModel;

typedef struct {
  int nstates;
  int ninputs;
  sfloat dt;
  Matrix* A;
  Matrix* B;
  Matrix* f;
  Matrix x0;
  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix);
  // int data_size;
} tiny_LtvModel;

void tiny_InitLtiModel(tiny_LtiModel* model);
void tiny_InitLtvModel(tiny_LtvModel* model);

typedef struct {
  Matrix x;   ///< state vector
  Matrix u;   ///< control input vector
  sfloat t;   ///< time
  sfloat dt;  ///< time step
} tiny_KnotPoint;

void tiny_InitKnotPoint(tiny_KnotPoint* z);

typedef struct {
  sfloat reg;
  sfloat reg_min;
  sfloat reg_max;
  sfloat penalty;
  sfloat penalty_max;
  sfloat penalty_mul;
  int max_outer_iters;
  int max_search_iters;
  sfloat riccati_tol;
  sfloat cstr_tol;
} tiny_Solver;

void tiny_InitSolver(tiny_Solver* solver);

typedef struct tiny_ProblemData {
  int nstates;
  int ninputs;
  int nhorizon;
  int ncstr_states;
  int ncstr_inputs;
  int ncstr_goal;
  Matrix Q;
  Matrix R;
  Matrix q;
  Matrix r;
  Matrix Qf;
  Matrix qf;
  Matrix* X_ref;
  Matrix* U_ref;
  sfloat dt;
  Matrix x0;
  Matrix* K;
  Matrix* d;
  Matrix* P;
  Matrix* p;
  Matrix* input_duals;
  Matrix* state_duals;
  Matrix goal_dual;
  Matrix Acstr_state;
  Matrix bcstr_state;
  Matrix Acstr_input;
  Matrix bcstr_input;
  // int data_size;
} tiny_ProblemData;

void tiny_InitProblemData(tiny_ProblemData* prob);