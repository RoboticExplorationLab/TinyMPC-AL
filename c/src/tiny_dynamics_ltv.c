#include "tiny_dynamics_ltv.h"

void tiny_DynamicsLtv(Matrix* xn, const Matrix x, const Matrix u,
                      const tiny_LtvModel model, const int k) {
  slap_MatrixCopy(*xn, model.f[k]);                        
  slap_MatMulAdd(*xn, model.A[k], x, 1, 1);      // x[k+1] += A * x[k]
  slap_MatMulAdd(*xn, model.B[k], u, 1, 1);      // x[k+1] += B * u[k]
}

// Roll out the closed-loop dynamics with K and d from backward pass,
// provided function to calculate Ak, Bk. Calculate new X, U in place
enum slap_ErrorCode tiny_ForwardPassLtv(
    Matrix* X, Matrix* U, const tiny_ProblemData prob, 
    const tiny_LtvModel model) {
  int N = prob.nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    // Control input: u = - d - K*x
    slap_MatrixCopy(U[k], prob.d[k]);              // u[k] = -d[k]
    slap_MatMulAdd(U[k], prob.K[k], X[k], -1, -1);   // u[k] -= K[k] * x[k]
    // Next state: x = A*x + B*u + f
    tiny_DynamicsLtv(&X[k + 1], X[k], U[k], model, k);
  }
  return SLAP_NO_ERROR;
}

void tiny_UpdateHorizonJacobians(tiny_LtvModel* model, tiny_ProblemData prob) {
    for (int i = 0; i < prob.nhorizon-1; ++i) {  
    model->get_jacobians(&(model->A[i]), &(model->B[i]), prob.X_ref[i], prob.U_ref[i]);
    tiny_Bicycle5dNonlinearDynamics(&(model->f[i]), prob.X_ref[i], prob.U_ref[i]);
    slap_MatMulAdd(model->f[i], model->A[i], prob.X_ref[i], -1, 1);
    slap_MatMulAdd(model->f[i], model->B[i], prob.U_ref[i], -1, 1);
  }   
}