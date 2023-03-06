#include "augmented_lagrangian_lqr.h"

const tiny_LinearDiscreteModel kDefaultLinearDiscreteModel = {
  .nstates = 0,
  .ninputs = 0,
  .dt = 0.0,
  .A = kNullMat,
  .B = kNullMat,
  .f = kNullMat,
  .x0 = kNullMat,
};

const tiny_KnotPoint kDefaultKnotPoint = {
  .x = kNullMat,
  .u = kNullMat,
  .t = 0.0,
  .dt = 0.0,
};

const tiny_Solver kDefaultSolver = {
  .regu = 0.0,
  .input_duals = NULL,
  .state_duals = NULL,
  .goal_dual = kNullMat,
  .penalty_min = 0.0,
  .penalty_max = 0.0,
  .penalty_mul = 0.0,
};

const tiny_ProblemData kDefaultProblemData = {
  .nstates = 0,
  .ninputs = 0,
  .nhorizon = 0,
  .ncstr_states = 0,
  .ncstr_inputs = 0,
  .ncstr_goal = 0,
  .Q = kNullMat,
  .R = kNullMat,
  .q = kNullMat,
  .r = kNullMat,
  .Qf = kNullMat,
  .u_max = kNullMat,
  .u_min = kNullMat,
  .x_max = kNullMat,
  .x_min = kNullMat,
  .X_ref = NULL,
  .U_ref = NULL,
  .dt = 0.0,
  .x0 = kNullMat,
  .K = NULL,
  .d = NULL,
  .P = NULL,
  .p = NULL,  
};

enum slap_ErrorCode tiny_BackwardPass(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    const Matrix* X, const Matrix* U, const tiny_Solver solver, Matrix S_temp)
{
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_ForwardPass(
    const tiny_ProblemData prob, const tiny_LinearDiscreteModel model, 
    Matrix* X, Matrix* U, const tiny_Solver solver)
{
  return SLAP_NO_ERROR;
}

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
    tiny_ProblemData prob, const tiny_LinearDiscreteModel model,
    Matrix* X, Matrix* U, const int verbose)
{
  return SLAP_NO_ERROR;
}

// [x-p.x_max;-x + p.x_min]
Matrix tiny_IneqInputs(
    const tiny_ProblemData prob, const Matrix u) 
{
  double ineq_data[prob.ncstr_inputs];
  Matrix mat = slap_MatrixFromArray(prob.ncstr_inputs, 1, ineq_data);
  Matrix upper_half = slap_CreateSubMatrix(mat, 0, 0, prob.ninputs, 0);
  Matrix lower_half = slap_CreateSubMatrix(mat, prob.ninputs, 0, 
                                          prob.ncstr_inputs, 0);
  slap_MatrixAddition(lower_half, u, prob.u_max, -1);
  return kNullMat;
}

Matrix tiny_IneqInputsJacobian(
    const tiny_ProblemData prob, const Matrix u) 
{
  return kNullMat;
}
  
Matrix tiny_ActiveIneqMask(
    const tiny_Solver solver, const Matrix ineq) 
{    
  return kNullMat;
}