#include "utils.h"

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
  Matrix* K;
  Matrix* d;
  Matrix* P;
  Matrix* p;
  Matrix* input_duals;
  Matrix* state_duals;
  Matrix goal_dual;
  int data_size;
} tiny_ProblemData;

enum tiny_ErrorCode tiny_SetProblemDims(tiny_ProblemData* prob, const int nstates, 
    const int ninputs, const int nhorizon);

enum tiny_ErrorCode tiny_AddInputLimitConstraints(tiny_ProblemData* prob, 
    sfloat* umin, sfloat* umax);

enum tiny_ErrorCode tiny_AddStateLimitConstraints(tiny_ProblemData* prob, 
    sfloat* xmin, sfloat* xmax);

