#include "data_struct.h"

void tiny_InitProblemData(tiny_ProblemData* prob);

void tiny_AddStageCost(sfloat* cost, const tiny_ProblemData prob,
                       const Matrix x, const Matrix u, const int k);

void tiny_AddTerminalCost(sfloat* cost, const tiny_ProblemData prob,
                          const Matrix x);

void tiny_ExpandStageCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                          Matrix* hes_el_uu, Matrix* grad_el_u,
                          const tiny_ProblemData prob, const Matrix x,
                          const Matrix u, const int k);

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob, const Matrix x);

enum slap_ErrorCode tiny_BackwardPassLti(tiny_ProblemData* prob,
                                         const tiny_LtiModel model,
                                         const tiny_Solver solver,
                                         const Matrix* X, const Matrix* U,
                                         Matrix G_temp);

enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_LtiModel model, const tiny_Solver solver,
    const Matrix* X, const Matrix* U, Matrix* G_temp, Matrix* ineq_temp);

enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtiModel model);

enum slap_ErrorCode tiny_AugmentedLagrangianLqr(Matrix* X, Matrix* U,
                                                tiny_ProblemData* prob,
                                                tiny_Solver* solver,
                                                const tiny_LtiModel model,
                                                const int verbose);

sfloat tiny_RiccatiConvergence(const tiny_ProblemData prob);

void tiny_IneqInputs(Matrix* ineq, const tiny_ProblemData prob, const Matrix u);

void tiny_IneqInputsJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

void tiny_IneqStates(Matrix* ineq_state, const tiny_ProblemData prob,
                     const Matrix x);

void tiny_IneqStatesJacobian(Matrix* ineq_jac, const tiny_ProblemData prob);

void tiny_ActiveIneqMask(Matrix* mask, const Matrix input_dual,
                         const Matrix ineq);

void tiny_ClampIneqDuals(Matrix* dual, const Matrix new_dual);

void tiny_DynamicsLti(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtiModel model);