#pragma once

#ifdef __cplusplus
 extern "C" {
#endif

#include <math.h>

#include <slap_arduino.h>

#define Hdyn 0.1 // dt for dynamics

/*
state: position x, y
input: velocity xdot, ydot
forward Euler integration
*/

void tiny_PointMassDynamics_Raw(sfloat* xn, const sfloat* x,
                                         const sfloat* u);

void tiny_PointMassDynamics(Matrix* xn, const Matrix x,
                                     const Matrix u);

void tiny_PointMassGetJacobianA_Raw(sfloat* A, const sfloat* x,
                                    const sfloat* u);

void tiny_PointMassGetJacobianB_Raw(sfloat* B, const sfloat* x,
                                    const sfloat* u);

void tiny_PointMassGetJacobians(Matrix* A, Matrix* B, const Matrix x,
                                const Matrix u);
                                
#ifdef __cplusplus
}
#endif                                
