#include "augmented_lagrangian_lqr.h"

enum slap_ErrorCode tiny_BackwardPass(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    const tiny_KnotPoint* Z, const tiny_Solver solver, Matrix S_temp)
{
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_ForwardPass(
    const tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    tiny_KnotPoint* Z, const tiny_Solver solver)
{
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model,
    tiny_KnotPoint* Z, const int verbose)
{
  return SLAP_NO_ERROR;
}

Matrix slap_IneqInputs(
    const tiny_ProblemData prob, const tiny_KnotPoint* Z) 
{
  return kNullMat;
}

Matrix slap_IneqInputsJacobian(
    const tiny_ProblemData prob, const tiny_KnotPoint* Z) 
{
  return kNullMat;
}
  
Matrix slap_ActiveIneqMask(
    const tiny_Solver solver, const Matrix ineq) 
{    
  return kNullMat;
}