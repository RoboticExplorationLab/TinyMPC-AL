#pragma once

#include <math.h>

#include "slap/slap.h"

//========================================
// Bicycle model parameters
// X = [x; y; theta; v; delta] : x, y, yaw, linear vel, steering angle
// U = [a; delta_dot] : linear accel and steering rate
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
void tiny_Bicycle5dNonlinearDynamics_Raw(double* xn, const double* x,
                                         const double* u);

void tiny_Bicycle5dNonlinearDynamics(Matrix* xn, const Matrix x,
                                     const Matrix u);

//========================================
// Codes generated from julia/bicycle_tvlqr
// Jacobians of discrete dynamics of bicycle model with predefined model params
//========================================
void tiny_Bicycle5dGetJacobianA_Raw(double* A, const double* x,
                                    const double* u);

void tiny_Bicycle5dGetJacobianB_Raw(double* B, const double* x,
                                    const double* u);

void tiny_Bicycle5dGetJacobians(Matrix* A, Matrix* B, const Matrix x,
                                const Matrix u);