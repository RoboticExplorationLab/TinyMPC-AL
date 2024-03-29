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


// TODO (sschoedel): replace tiny_LtiModel and tiny_LtvModel with tiny_Model
typedef struct {
  int nstates;
  int ninputs;
  double dt;
  Matrix A;
  Matrix B;
  Matrix f;
  Matrix x0;
  // int data_size;  ///< number of doubles need to store the data
} tiny_LtiModel;

typedef struct {
  int nstates;
  int ninputs;
  double dt;
  Matrix* A;
  Matrix* B;
  Matrix* f;
  Matrix x0;
  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix);
  // int data_size;
} tiny_LtvModel;

enum modelType {kTimeInvariant = 0, kTimeVariant = 1};

typedef struct {
  enum modelType type;
  int nstates;
  int ninputs;
  double dt;
  Matrix* A;
  Matrix* B;
  Matrix* f;
  Matrix x0;
  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonlinear_dynamics)(Matrix*, const Matrix, const Matrix);
  // int data_size;
} tiny_Model;

void tiny_InitLtiModel(tiny_LtiModel* model);
void tiny_InitLtvModel(tiny_LtvModel* model);

typedef struct {
  Matrix x;   ///< state vector
  Matrix u;   ///< control input vector
  double t;   ///< time
  double dt;  ///< time step
} tiny_KnotPoint;

void tiny_InitKnotPoint(tiny_KnotPoint* z);

typedef struct {
  double reg;
  double reg_min;
  double reg_max;
  double penalty;
  double penalty_max;
  double penalty_mul;
  int max_primal_iters;
  int max_search_iters;
  double riccati_tol;
  double cstr_tol;
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
  Matrix u_max;
  Matrix u_min;
  Matrix x_max;
  Matrix x_min;
  Matrix* X_ref;
  Matrix* U_ref;
  double dt;
  Matrix x0;
  Matrix* K;
  Matrix* d;
  Matrix* P;
  Matrix* p;
  Matrix* input_duals;
  Matrix* state_duals;
  Matrix goal_dual;
  // int data_size;
} tiny_ProblemData;

void tiny_InitProblemData(tiny_ProblemData* prob);