#include "data_struct.h"

void tiny_DynamicsLtv(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtvModel model, const int k);

enum slap_ErrorCode tiny_ForwardPassLtv(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtvModel model);

void tiny_UpdateHorizonJacobians(tiny_LtvModel* model, tiny_ProblemData prob);