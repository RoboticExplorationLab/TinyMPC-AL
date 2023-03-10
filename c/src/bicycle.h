#pragma once

#include <math.h>

#include "slap/slap.h"

//========================================
// Bicycle model parameters
//========================================
// struct tiny_Model_Bicycle {
//   double drive_min[2];
//   double drive_max[2];
//   double u_min[2];
//   double u_max[2];
// } tiny_DefaultModel_Bicycle = {{-2, -0.5}, {2, 0.5}, {-4, -0.7}, {4, 0.7}};

//========================================
// Codes generated from julia/bicycle_tvlqr
// Discrete dynamics of bicycle model with predefined model params
//========================================
void tiny_DiscreteDynamics_Raw(double* xn, const double* x, const double* u);

void tiny_DiscreteDynamics(Matrix* xn, const Matrix x, const Matrix u);

//========================================
// Codes generated from julia/bicycle_tvlqr
// Jacobians of discrete dynamics of bicycle model with predefined model params
//========================================
void tiny_GetJacobianA_Raw(double* A, const double* x, const double* u);

void tiny_GetJacobianB_Raw(double* B, const double* x, const double* u);

void tiny_GetJacobians(Matrix* A, Matrix* B, const Matrix x, const Matrix u);