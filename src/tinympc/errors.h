#pragma once

#include "slap/slap.h"

// Inherit from slap_ErrorCode and expand new errors for tinympc
enum tiny_ErrorCode {
  TINY_SLAP_ERROR = 0,
  TINY_MATRIX_NOT_PSD,
  TINY_PROBLEM_INFEASIBLE,
  TINY_NO_ERROR,
};
