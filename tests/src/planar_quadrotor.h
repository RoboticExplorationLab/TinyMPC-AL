#pragma once

#include <math.h>

#include "slap/slap.h"

//========================================
// Planar quadrotor model parameters
//========================================
// struct tiny_Model_PlanarQuadrotor {
//   sfloat g;
//   sfloat m;
//   sfloat l;
//   sfloat J;
//   sfloat umin[2];
//   sfloat umax[2];
// } tiny_DefaultModel_PlanarQuadrotor = {9.81, 1,      0.018,
//                                        0.3,  {0, 0}, {19.62, 19.62}};

//========================================
// Codes generated from julia/planar_quad_gen
// Discrete dynamics of planar quadrotor
//========================================
void tiny_PQuadNonlinearDynamics_Raw(sfloat* xn, const sfloat* x,
                                     const sfloat* u);

void tiny_PQuadNonlinearDynamics(Matrix* xn, const Matrix x, const Matrix u);

//========================================
// Codes generated from julia/planar_quad_gen
// Jacobians of discrete dynamics of planar quadrotor
//========================================
void tiny_PQuadGetJacobianA_Raw(sfloat* A, const sfloat* x, const sfloat* u);

void tiny_PQuadGetJacobianB_Raw(sfloat* B, const sfloat* x, const sfloat* u);

void tiny_PQuadGetJacobians(Matrix* A, Matrix* B, const Matrix x,
                            const Matrix u);