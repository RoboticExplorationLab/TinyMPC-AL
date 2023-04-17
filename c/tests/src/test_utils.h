#pragma once
#include "slap/slap.h"
// #include "riccati/riccati_solver.h"
sfloat SumOfSquaredError(const sfloat* x, const sfloat* y, int len);

sfloat SumOfSquaredErrorMatrices(const sfloat* x, Matrix* Y, const int num);

// void DiscretesfloatIntegratorDynamics(sfloat h, sfloat dim, Matrix* A,
// Matrix* B);

// RiccatiSolver* sfloatIntegratorProblem();
