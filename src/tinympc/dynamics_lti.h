#include "slap/slap.h"
#include "data_struct.h"

void tiny_DynamicsLti(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtiModel model);

enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtiModel model);