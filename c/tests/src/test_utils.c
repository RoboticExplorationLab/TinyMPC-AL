#include "test_utils.h"

#include <math.h>

// #include "riccati/riccati_solver.h"
#include "slap/matrix.h"

sfloat SumOfSquaredError(const sfloat* x, const sfloat* y, int len) {
  sfloat err = 0;
  for (int i = 0; i < len; ++i) {
    sfloat diff = x[i] - y[i];
    err += diff * diff;
  }
  return sqrt(err);
}

sfloat SumOfSquaredErrorMatrices(const sfloat* x, Matrix* Y, const int num) {
  sfloat err = 0;
  int k = 0;
  for (int i = 0; i < num; ++i) {
    for (int j = 0; j < Y[i].cols * Y[i].rows; ++j) {
      sfloat diff = x[k++] - Y[i].data[j];
      err += diff * diff;
    }
  }
  return sqrt(err);
}

// void DiscretesfloatIntegratorDynamics(sfloat h, sfloat dim, Matrix* A,
// Matrix* B) {
//   int nstates = 2 * dim;
//   slap_MatrixSetConst(A, 0.0);
//   slap_MatrixSetConst(B, 0.0);
//   for (int i = 0; i < nstates; ++i) {
//     slap_MatrixSetElement(A, i, i, 1.0);
//   }
//   sfloat b = h * h / 2;
//   for (int i = 0; i < dim; ++i) {
//     slap_MatrixSetElement(A, i, i + dim, h);
//     slap_MatrixSetElement(B, i, i, b);
//     slap_MatrixSetElement(B, i + dim, i, h);
//   }
// }

// RiccatiSolver* sfloatIntegratorProblem() {
//   const int dim = 2;
//   const int nstates = 2 * dim;
//   const int ninputs = dim;
//   const int nhorizon = 11;
//   const sfloat h = 0.1;

//   // Generate the problem
//   RiccatiSolver* solver = ulqr_NewRiccatiSolver(nstates, ninputs, nhorizon);

//   // Generate dynamics matrices and copy into problem
//   Matrix A = slap_NewMatrix(nstates, nstates);
//   Matrix B = slap_NewMatrix(nstates, ninputs);
//   DiscretesfloatIntegratorDynamics(h, dim, &A, &B);
//   const sfloat f[4] = {
//       1,
//       -1,
//       0,
//       0,
//   };
//   ulqr_SetDynamics(solver, A.data, B.data, f, 0, nhorizon - 1);

//   return solver;
// }
