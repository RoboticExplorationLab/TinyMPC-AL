#pragma once

#include "slap/slap.h"
#include "data_struct.h"

void tiny_AddStageCost(sfloat* cost, const tiny_ProblemData prob,
                       const Matrix x, const Matrix u, const int k);

void tiny_AddTerminalCost(sfloat* cost, const tiny_ProblemData prob,
                          const Matrix x);

void tiny_ExpandStageCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                          Matrix* hes_el_uu, Matrix* grad_el_u,
                          const tiny_ProblemData prob, const int k);

void tiny_ExpandTerminalCost(Matrix* hes_el_xx, Matrix* grad_el_x,
                             const tiny_ProblemData prob);