#pragma once

#ifdef __cplusplus
extern "C" {
#endif


// #include "riccati/riccati_solver.h"
double SumOfSquaredError(const double* x, const double* y, int len);

// void DiscreteDoubleIntegratorDynamics(double h, double dim, Matrix* A,
// Matrix* B);

// RiccatiSolver* DoubleIntegratorProblem();

#ifdef __cplusplus
}
#endif