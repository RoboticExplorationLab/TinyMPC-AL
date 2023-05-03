#ifndef TYPES_H
# define TYPES_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "slap/slap.h"
#include "constants.h"
#include "utils.h"
#include "errors.h"

// for a horizon of N x(0)->x(N-1), need N-1 matrices
typedef struct {
  int nstates;
  int ninputs;
  int nhorizon;

  int ltv;            ///< Boolean, true if model is LTV  
  int affine;         ///< Boolean, true if model is affine
  sfloat dt;          ///< Sample time Ts of the discrete model

  Matrix* A;
  Matrix* B;
  Matrix* f;

  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonl_model)(Matrix*, const Matrix, const Matrix);
  int data_size;
} tiny_Model;

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
  int data_size;
} tiny_Solution;


/**
 * Solver return information
 */
typedef struct {
  int iter_al;        ///< Number of AL iterations taken
  int iter_riccati;   ///< Number of Riccati iterations taken
  int status_val;     ///< Integer, status defined in constants.h

  sfloat obj_pri;     ///< primal objective
  sfloat obj_al;      ///< Augmented Lagrangian objective
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
  sfloat reg_mul;             ///< Regularization update multiplier
  int    en_reg_update;       ///< Boolean, enable regularization update (tighter solve)
  
  sfloat penalty_init;        ///< Initial penalty
  sfloat penalty_max;         ///< Maximum penalty
  sfloat penalty_mul;         ///< Penalty multiplier

  sfloat alpha_mul;           ///< Line-search step multiplier

  int    max_iter_al;         ///< Maximum number of AL iterations
  int    max_iter_riccati;    ///< Maximum number of Riccati solve iterations
  int    max_iter_ls;         ///< Maximum number of line-search iterations

  sfloat tol_abs_riccati;     ///< Riccati solve tolerance
  sfloat tol_abs_cstr;        ///< Constraint tolerance

  int    en_cstr_states;      ///< Boolean, enable inequality constraints on states
  int    en_cstr_inputs;      ///< Boolean, enable inequality constraints on inputs
  int    en_cstr_goal;        ///< Boolean, enable equality constraint on goal

  int    verbose;             ///< Integer, level to write out progress
  int    adaptive_horizon;    ///< Integer, after `adaptive_horizon` steps, use the second model with longer interval; if 0, disabled 
  int    check_riccati;       ///< Boolean, if 0, then termination checking is disabled
  int    check_al;            ///< Boolean, if 0, then termination checking is disabled
  int    warm_start;          ///< boolean, enable warm start
  sfloat time_limit;          ///< Time limit of each MPC step; if 0, disabled
} tiny_Settings;

// void tiny_InitSettings(tiny_Settings* solver);


/**
 * Data structure
 */
typedef struct {
  tiny_Model* model;    ///< System model
  Matrix x0;

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
  int data_size;
} tiny_Data;

// void tiny_InitProblemData(tiny_ProblemData* prob);

typedef struct {
  tiny_Data*        data;      ///< problem data
  tiny_Settings*    stgs;      ///< problem settings
  tiny_Solution*    soln;      ///< problem solution
  tiny_Info*        info;      ///< solver information

  sfloat reg;
  sfloat alpha;
  sfloat penalty;
  // Temporary data
  Matrix Q_temp;
  Matrix c_temp;

  Matrix Qxx;
  Matrix Qxu;
  Matrix Qux;
  Matrix Quu;
  Matrix Qx;
  Matrix Qu;
  
  Matrix cu;
  Matrix cu2;
  Matrix cu_jac;
  Matrix cu_jac2;
  Matrix cu_mask;
  Matrix YU_hat;

  Matrix cx;
  Matrix cx2;
  Matrix cx_jac;
  Matrix cx_jac2;
  Matrix cx_mask;
  Matrix YX_hat;

  Matrix cg;

  int data_size;      ///< sum data size of all temporary data //TODO: + model + solution 
  int first_run;      ///< flag indicating whether the solve function has been run before
} tiny_Workspace;


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef TYPES_H
