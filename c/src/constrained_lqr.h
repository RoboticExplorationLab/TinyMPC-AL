#include "tiny_struct.h"
#include "tiny_cost.h"
#include "tiny_dynamics.h"

enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_Solver solver, 
    const tiny_LtiModel model, const Matrix* X, const Matrix* U, 
    Matrix* Q_temp, Matrix* ineq_temp);
    
enum slap_ErrorCode tiny_MpcLti(
    Matrix* X, Matrix* U, tiny_ProblemData* prob, tiny_Solver* solver,
    const tiny_LtiModel model, const int verbose);

enum slap_ErrorCode tiny_MpcLtv(
    Matrix* X, Matrix* U, tiny_ProblemData* prob, tiny_Solver* solver,
    tiny_LtiModel* model, const int verbose, 
    void (*get_jacobians)(Matrix, Matrix, const Matrix, const Matrix));

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