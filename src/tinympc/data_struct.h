#pragma once

#include "utils.h"

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