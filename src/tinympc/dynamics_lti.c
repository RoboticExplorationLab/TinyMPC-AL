#include "dynamics_lti.h"

void tiny_DynamicsLti(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtiModel model) {
  slap_Copy(*xn, model.f);
  slap_MatMulAdd(*xn, model.A, x, 1, 1);  // x[k+1] += A * x[k]
  slap_MatMulAdd(*xn, model.B, u, 1, 1);  // x[k+1] += B * u[k]
}

// Roll out the closed-loop dynamics with K and d from backward pass,
// calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPassLti(Matrix* X, Matrix* U,
                                        const tiny_ProblemData prob,
                                        const tiny_LtiModel model) {
  int N = prob.nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    // delta_x and delta_u over previous X, U
    // Control input: u = uf - d - K*(x - xf)
    slap_Copy(U[k], prob.d[k]);               // u[k] = d[k]
    slap_MatMulAdd(U[k], prob.K[k], X[k], -1, -1);  // u[k] += K[k] * x[k]
    // Next state: x = A*x + B*u + f
    tiny_DynamicsLti(&X[k + 1], X[k], U[k], model);
  }
  return SLAP_NO_ERROR;
}