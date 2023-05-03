#ifndef CONSTRAINT_LINEAR_H
# define CONSTRAINT_LINEAR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"

enum tiny_ErrorCode tiny_EvalInputConstraint(tiny_Workspace* work, const int k);

// void tiny_EvalInputConstraintOffset(Matrix* cu, const tiny_ProblemData prob);

// void tiny_EvalInputConstraintJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

enum tiny_ErrorCode tiny_EvalStateConstraint(tiny_Workspace* work, const int k);

// void tiny_EvalStateConstraintOffset(Matrix* cx, const tiny_ProblemData prob);

// void tiny_EvalStateConstraintJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

enum tiny_ErrorCode tiny_ActiveIneqMask(Matrix* mask, const Matrix dual,
                         const Matrix eval);

enum tiny_ErrorCode tiny_ProjectOrthantDuals(Matrix* dual, const Matrix new_dual);

int IsConstrained(tiny_Workspace* work);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CONSTRAINT_LINEAR_H