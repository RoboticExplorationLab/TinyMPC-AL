#include "tiny_struct.h"
#include "tiny_cost.h"
#include "tiny_constraint.h"
#include "tiny_dynamics_lti.h"

enum slap_ErrorCode tiny_ConstrainedBackwardPassLti(
    tiny_ProblemData* prob, const tiny_Solver solver, 
    const tiny_LtiModel model, const Matrix* X, const Matrix* U, 
    Matrix* Q_temp, Matrix* ineq_temp);
    
enum slap_ErrorCode tiny_MpcLti(
    Matrix* X, Matrix* U, tiny_ProblemData* prob, tiny_Solver* solver,
    const tiny_LtiModel model, const int verbose);

