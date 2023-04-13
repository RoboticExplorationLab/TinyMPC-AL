#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "slap/slap.h"

sfloat SumOfSquaredError(const sfloat* x, const sfloat* y, const int len);

sfloat SumOfSquaredErrorMatrices(const sfloat* x, Matrix* Y, const int num);

#ifdef __cplusplus
}
#endif