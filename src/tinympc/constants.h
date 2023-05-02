#ifndef CONSTANTS_H
# define CONSTANTS_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


/******************
* Solver Status  *
******************/
# define TINY_SOLVED (1)
# define TINY_MAX_ITER_AL_REACHED (-2)
# define TINY_MAX_ITER_RICCATI_REACHED (-3)
# define TINY_MAX_ITER_LS_REACHED (-4)
# define TINY_UNSOLVED (-10)            /* Unsolved. Only setup function has been called */


/******************
* Solver Errors  *
******************/
// Inherit from slap_ErrorCode and expand new errors for tinympc
enum tiny_ErrorCode {
  TINY_SLAP_ERROR = 0,
  TINY_MATRIX_NOT_PSD,
  TINY_PROBLEM_INFEASIBLE,
  TINY_NO_ERROR,
};


/**********************************
* Solver Parameters and Settings *
**********************************/

# define REG_MIN (1e-6)
# define REG_MAX (1e2)
# define REG_MUL (1.6)
# define EN_REG_UPDATE (0)

# define PENALTY_INIT (1e0)
# define PENALTY_MAX (1e6)
# define PENALTY_MUL (10.0)

# define ALPHA_MUL (0.5)

# define MAX_ITER_AL (10)
# define MAX_ITER_RICCATI (10)
# define MAX_ITER_LS (10)

# define TOL_ABS_RICCATI (1e-2)
# define TOL_ABS_CSTR (1e-2)

# define EN_CSTR_STATES (1)
# define EN_CSTR_INPUTS (1)
# define EN_CSTR_GOAL (1)

# define VERBOSE (1)
# define ADAPTIVE_HORIZON (0)
# define CHECK_TERMINATION (0)
# define WARM_START (1)
# define TIME_LIMIT (0.0)


# ifndef TINY_NULL_MAT
#  define TINY_NULL_MAT  \
  ((Matrix){      \
      0,          \
      0,          \
      0,          \
      0,          \
      TINY_NULL,       \
      slap_DENSE, \
  })
# endif /* ifndef TINY_NULL_MAT */

# ifndef TINY_NULL
#  define TINY_NULL 0
# endif /* ifndef TINY_NULL */

# ifndef TINY_NAN
#  define TINY_NAN ((sfloat)0x7fc00000UL)  // not a number
# endif /* ifndef TINY_NAN */

# ifndef TINY_INFTY
#  define TINY_INFTY ((sfloat)1e30)        // infinity
# endif /* ifndef TINY_INFTY */


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CONSTANTS_H
