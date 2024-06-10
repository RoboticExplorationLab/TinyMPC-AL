#include "pointmass.h"

/*
state: position x, y
input: velocity xdot, ydot
forward Euler integration
*/

void tiny_PointMassDynamics_Raw(sfloat* xn, const sfloat* x, const sfloat* u) {
  xn[0] = x[0] + Hdyn * u[0];
  xn[1] = x[1] + Hdyn * u[1];
}

void tiny_PointMassDynamics(Matrix* xn, const Matrix x,
                                     const Matrix u) {
  tiny_PointMassDynamics_Raw(xn->data, x.data, u.data);
}

void tiny_PointMassGetJacobianA_Raw(sfloat* A, const sfloat* x, const sfloat* u) {
  A[0] = 1;
  A[1] = 0;
  A[2] = 0;
  A[3] = 1;
}

void tiny_PointMassGetJacobianB_Raw(sfloat* B, const sfloat* x, const sfloat* u) {
  B[0] = Hdyn;
  B[1] = 0;
  B[2] = 0;
  B[3] = Hdyn;
}

void tiny_PointMassGetJacobians(Matrix* A, Matrix* B, const Matrix x,
                                const Matrix u) {
  tiny_PointMassGetJacobianA_Raw(A->data, x.data, u.data);
  tiny_PointMassGetJacobianB_Raw(B->data, x.data, u.data);
}