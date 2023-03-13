#pragma once

#include "slap/slap.h"

enum slap_ErrorCode get_augmented_lagrangian_cost(const Matrix A, const Matrix b,
                            const Matrix Q, const Matrix R,
                            Matrix *x, Matrix *u, int *L, Matrix *L_grad, Matrix *L_hess,
                            const int NSTATES, const int NINPUTS, const int NHORIZON,
                            Matrix xf, Matrix q, Matrix r, Matrix *Khist,
                            Matrix *dhist, Matrix *Phist, Matrix *phist, Matrix *xhist,
                            Matrix *uhist, Matrix yhist, Matrix x0, const int rho);