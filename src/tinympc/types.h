#ifndef TYPES_H
# define TYPES_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "slap/slap.h"


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
  int nhorizon;
  sfloat dt;
  Matrix* A;
  Matrix* B;
  Matrix* f;
  Matrix x0;
  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonl_model)(Matrix*, const Matrix, const Matrix);
  // int data_size;
} tiny_LtvModel;

void tiny_InitLtiModel(tiny_LtiModel* model);
void tiny_InitLtvModel(tiny_LtvModel* model);


/**
 * Solution structure
 */
typedef struct {
  Matrix* X;      ///< State trajectory solution 
  Matrix* U;      ///< Input trajectory solution

  Matrix* YX;     ///< Input duals associated to \f$l <= A*X[k] <= u\f$
  Matrix* YU;     ///< State duals associated to \f$l <= A*U[k] <= u\f$
  Matrix  YG;     ///< Goal dual associated to \f$X[N] = xg\f$

  Matrix* K;
  Matrix* d;
  Matrix* P;
  Matrix* p;
} tiny_Solution;


/**
 * Solver return information
 */
typedef struct {
  int iter;           ///< Number of iterations taken
  int status_val;     ///< Integer, status defined in constants.h

  sfloat obj_val;     ///< primal objective
  sfloat pri_res;     ///< norm of primal residual
  sfloat dua_res;     ///< norm of dual residual
} tiny_Info;


/**********************************
* Main structures and Data Types *
**********************************/

/**
 * Settings struct
 */
typedef struct {
  sfloat reg_min;             ///< Minimum regularization
  sfloat reg_max;             ///< Maximum regularization
  int    en_reg_update;       ///< Enable regularization update (tighter solve)
  
  sfloat penalty_init;        ///< Initial penalty
  sfloat penalty_max;         ///< Maximum penalty
  sfloat penalty_mul;         ///< Penalty multiplier

  sfloat alpha_mul;           ///< Line-search step multiplier

  int    max_iter_al;         ///< Maximum number of AL iterations
  int    max_iter_riccati;    ///< Maximum number of Riccati solve iterations
  int    max_iter_ls;         ///< Maximum number of line-search iterations

  sfloat tol_abs_riccati;     ///< Riccati solve tolerance
  sfloat tol_abs_cstr;        ///< Constraint tolerance

  int    en_cstr_states;      ///< Enable inequality constraints on states
  int    en_cstr_inputs;      ///< Enable inequality constraints on inputs
  int    en_cstr_goal;        ///< Enable equality constraint on goal

  int    adaptive_horizon;    ///< Integer, from this step, use the second model with longer step; if 0, disabled 
  
  sfloat time_limit;          ///< Time limit of each MPC step; if 0, disabled
} tiny_Settings;

void tiny_InitSettings(tiny_Settings* solver);


/**
 * Data structure
 */
typedef struct {
  int nstates;
  int ninputs;
  int nhorizon;

  Matrix  Q;
  Matrix  R;
  Matrix  Qf;
  Matrix* q;
  Matrix* r;
  Matrix  qf;

  Matrix* X_ref;
  Matrix* U_ref;

  Matrix Acx;
  Matrix bcx;
  Matrix Acu;
  Matrix bcu;
} tiny_ProblemData;

void tiny_InitProblemData(tiny_ProblemData* prob);

typedef struct {
  tiny_ProblemData* data;      ///< problem data
  tiny_Settings*    stgs;      ///< problem settings
  tiny_Solution*    soln;      ///< problem solution
  tiny_Info*        info;      ///< solver information

  Matrix Qxx;
  Matrix Qxu;
  Matrix Qux;
  Matrix Quu;
  Matrix Qx;
  Matrix Qu;
  
  c_int first_run;      ///< flag indicating whether the solve function has been run before
} tiny_Workspace;


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef TYPES_H
