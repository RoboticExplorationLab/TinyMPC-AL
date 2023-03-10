#include "tiny_struct.h"

void tiny_InitProblemData(tiny_ProblemData* prob);

void tiny_AddStageCost(double* cost, const tiny_ProblemData prob,
                       const Matrix x, const Matrix u, const int k);

void tiny_AddTerminalCost(double* cost, const tiny_ProblemData prob,
                          const Matrix x);

void tiny_ExpandStageCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                          Matrix* hes_el_uu, Matrix* grad_el_u,
                          const tiny_ProblemData prob, const int k);

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob);

enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_Solver solver,
                                         const tiny_LinearDiscreteModel model,
                                         const Matrix* X, const Matrix* U,
                                         Matrix Q_temp);

// enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
//     tiny_ProblemData* prob, const tiny_Solver solver, 
//     const tiny_LinearDiscreteModel model, const Matrix* X, const Matrix* U, 
//     Matrix* Q_temp, Matrix* ineq_temp);

enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LinearDiscreteModel model);

// enum slap_ErrorCode tiny_AugmentedLagrangianLqr(
//     Matrix* X, Matrix* U, tiny_ProblemData* prob, tiny_Solver* solver,
//     const tiny_LinearDiscreteModel model, const int verbose);

double tiny_RiccatiConvergence(const tiny_ProblemData prob);

void tiny_IneqInputs(Matrix* ineq, const tiny_ProblemData prob, const Matrix u);

void tiny_IneqInputsOffset(Matrix* ineq_input, const tiny_ProblemData prob);

void tiny_IneqInputsJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

void tiny_IneqStates(Matrix* ineq_state, const tiny_ProblemData prob,
                     const Matrix x);

void tiny_IneqStatesOffset(Matrix* ineq_state, const tiny_ProblemData prob);

void tiny_IneqStatesJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

void tiny_ActiveIneqMask(Matrix* mask, const Matrix input_dual,
                         const Matrix ineq);

void tiny_ClampIneqDuals(Matrix* dual, const Matrix new_dual);

void tiny_DiscreteDynamics(Matrix* xn, const Matrix x, const Matrix u,
                           const tiny_LinearDiscreteModel model);