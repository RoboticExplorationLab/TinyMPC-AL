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
                                         const tiny_LtiModel model,
                                         Matrix Q_temp);

enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtiModel model);

enum slap_ErrorCode tiny_BackwardPassLtv(
    tiny_ProblemData* prob, const tiny_Solver solver,
    const tiny_LtvModel model, Matrix Q_temp);

enum slap_ErrorCode tiny_ForwardPassLtv(
    Matrix* X, Matrix* U, const tiny_ProblemData prob, 
    const tiny_LtvModel model);

// enum slap_ErrorCode tiny_BackwardPassLtvf(
//     tiny_ProblemData* prob, const tiny_Solver solver,
//     const tiny_LtvModel model, Matrix Q_temp);

// enum slap_ErrorCode tiny_ForwardPassLtvf(
//     Matrix* X, Matrix* U, const tiny_ProblemData prob, 
//     const tiny_LtvModel model);

void tiny_LtiDynamics(Matrix* xn, const Matrix x, const Matrix u,
                           const tiny_LtiModel model);

void tiny_LtvDynamics(Matrix* xn, const Matrix x, const Matrix u,
                           const tiny_LtvModel model, const int k);