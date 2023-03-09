#pragma once

#include "slap/slap.h"

#define kNullMat ((Matrix) {0, 0, 0, 0, NULL, slap_DENSE,})

typedef struct {
  int nstates;  
  int ninputs;
  double dt;
  Matrix A;
  Matrix B;
  Matrix f;
  Matrix x0;
} tiny_LinearDiscreteModel;

extern const tiny_LinearDiscreteModel kDefaultLinearDiscreteModel;

typedef struct {
  Matrix x;  ///< state vector
  Matrix u;  ///< control input vector
  double t;  ///< time
  double dt;  ///< time step
} tiny_KnotPoint;

extern const tiny_KnotPoint kDefaultKnotPoint; 

typedef struct {
  double regu;
  double regu_min;
  double regu_max;
  double penalty;
  double penalty_max;
  double penalty_mul;
  int max_primal_iters;
  int max_search_iters;
  double riccati_tol;
  double cstr_tol;
} tiny_Solver;

extern const tiny_Solver kDefaultSolver;

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
} tiny_ProblemData;

extern const tiny_ProblemData kDefaultProblemData;

void tiny_AddStageCost(
    double* cost, const tiny_ProblemData prob, 
    const Matrix x, const Matrix u, const int k);

void tiny_AddTerminalCost(
    double* cost, const tiny_ProblemData prob, const Matrix x);

void tiny_ExpandStageCost(
    Matrix* hes_el_xx, Matrix* grad_el_x, Matrix* hes_el_uu, Matrix* grad_el_u,
    const tiny_ProblemData prob, const Matrix x, const Matrix u, const int k);

void tiny_ExpandTerminalCost(
    Matrix* hes_el_xx, Matrix* grad_el_x, 
    const tiny_ProblemData prob, const Matrix x);

enum slap_ErrorCode tiny_BackwardPassLti(
    tiny_ProblemData* prob, const tiny_LinearDiscreteModel model, 
    const tiny_Solver solver, const Matrix* X, const Matrix* U,  Matrix G_temp);

enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_LinearDiscreteModel model, 
    const tiny_Solver solver, const Matrix* X, const Matrix* U,  
    Matrix* G_temp, Matrix* ineq_temp); 

enum slap_ErrorCode tiny_ForwardPassLti(
    Matrix* X, Matrix* U,
    const tiny_ProblemData prob, const tiny_LinearDiscreteModel model);

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
    Matrix* X, Matrix* U, tiny_ProblemData* prob, tiny_Solver* solver,
    const tiny_LinearDiscreteModel model, const int verbose); 

double tiny_RiccatiConvergence(tiny_ProblemData prob);

void tiny_IneqInputs(
  Matrix* ineq, const tiny_ProblemData prob, const Matrix u);

void tiny_IneqInputsJacobian(
    Matrix* ineq_jac, const tiny_ProblemData prob, const Matrix u);

void tiny_ActiveIneqMask(
  Matrix* mask, const Matrix input_dual, const Matrix ineq);    

void tiny_ClampIneqDuals(tiny_ProblemData* prob, Matrix new_dual, int k);

void tiny_DiscreteDynamics(
    Matrix* xn, const Matrix x, const Matrix u, 
    const tiny_LinearDiscreteModel model);