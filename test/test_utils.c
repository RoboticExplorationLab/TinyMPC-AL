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

// void DiscreteDoubleIntegratorDynamics(double h, double dim, Matrix* A,
// Matrix* B) {
//   int nstates = 2 * dim;
//   slap_MatrixSetConst(A, 0.0);
//   slap_MatrixSetConst(B, 0.0);
//   for (int i = 0; i < nstates; ++i) {
//     slap_MatrixSetElement(A, i, i, 1.0);
//   }
//   double b = h * h / 2;
//   for (int i = 0; i < dim; ++i) {
//     slap_MatrixSetElement(A, i, i + dim, h);
//     slap_MatrixSetElement(B, i, i, b);
//     slap_MatrixSetElement(B, i + dim, i, h);
//   }
// }

// RiccatiSolver* DoubleIntegratorProblem() {
//   const int dim = 2;
//   const int nstates = 2 * dim;
//   const int ninputs = dim;
//   const int nhorizon = 11;
//   const double h = 0.1;

//   // Generate the problem
//   RiccatiSolver* solver = ulqr_NewRiccatiSolver(nstates, ninputs, nhorizon);

//   // Generate dynamics matrices and copy into problem
//   Matrix A = slap_NewMatrix(nstates, nstates);
//   Matrix B = slap_NewMatrix(nstates, ninputs);
//   DiscreteDoubleIntegratorDynamics(h, dim, &A, &B);
//   const double f[4] = {
//       1,
//       -1,
//       0,
//       0,
//   };
//   ulqr_SetDynamics(solver, A.data, B.data, f, 0, nhorizon - 1);

//   return solver;
// }
