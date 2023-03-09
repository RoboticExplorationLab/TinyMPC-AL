/**
 * @file errors.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Utilities for the slap error handling system
 * @date 2022-01-30
 *
 * @copyright Copyright (c) 2022
 *
 * @addtogroup Utilities
 * @{
 */
#pragma once

enum slap_ErrorCode {
  SLAP_NO_ERROR = 0,
  SLAP_BAD_POINTER,
  SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  SLAP_BAD_MATRIX_DATA_POINTER,
  SLAP_INVALID_DIMENSION,
  SLAP_INVALID_STRIDE,
  SLAP_MATRIX_NOT_SQUARE,
  SLAP_MATRIX_NOT_DENSE,
  SLAP_CHOLESKY_FAIL,
  SLAP_INVALID_MATRIX,
  SLAP_INDEX_OUT_OF_BOUNDS,
};

const char* slap_ErrorString(enum slap_ErrorCode error_code);

#define SLAP_ERROR_STRING(error_code, ...)                                 \
  {                                                                        \
    fprintf(stderr, SLAP_COLOR_RED "slap Error %d: %s\n", (int)error_code, \
            slap_ErrorString(error_code));                                 \
    fprintf(stderr, SLAP_COLOR_RED "              ");                      \
    fprintf(stderr, SLAP_COLOR_RED __VA_ARGS__);                           \
    fprintf(stderr, SLAP_COLOR_RED " (%s:%d)\n", __FILE__, __LINE__);      \
  }

// TODO (brian): Add compile option to turn this on/off
#include <stdio.h>
#define SLAP_COLOR_RED "\x1b[31m"
#define SLAP_THROW_ERROR(error_code, message)                              \
  ((fprintf(stderr, SLAP_COLOR_RED "slap Error %d: %s\n", (int)error_code, \
            slap_ErrorString(error_code)),                                 \
    fprintf(stderr, SLAP_COLOR_RED "              %s (%s:%d)\n", message,  \
            __FILE__, __LINE__)),                                          \
   error_code)

#define SLAP_ASSERT(condition, error_code, return_value, ...) \
  if (!(condition)) {                                         \
    SLAP_ERROR_STRING(error_code, __VA_ARGS__)                \
    return return_value;                                      \
  }

#define SLAP_ASSERT_OK(error_code)                                     \
  SLAP_ASSERT((error_code) == SLAP_NO_ERROR, error_code, EXIT_FAILURE, \
              "Got Bad Return Code");

#define SLAP_ASSERT_VALID(mat, return_value, ...) \
  SLAP_ASSERT(slap_IsValid(mat), SLAP_INVALID_MATRIX, return_value, __VA_ARGS__)

#define SLAP_ASSERT_DENSE(mat, return_value, ...)                     \
  SLAP_ASSERT(slap_IsDense(mat), SLAP_MATRIX_NOT_DENSE, return_value, \
              __VA_ARGS__)

#define SLAP_ASSERT_SAME_SIZE(A, B, return_value, method_name)             \
  SLAP_ASSERT(                                                             \
      slap_NumRows(A) == slap_NumRows(B) &&                                \
          slap_NumCols(A) == slap_NumCols(B),                              \
      SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS, (return_value),                 \
      "%s: matrices must be the same size. Got sizes (%d,%d) and (%d,%d)", \
      method_name, slap_NumRows(A), slap_NumCols(A), slap_NumRows(B),      \
      slap_NumCols(B))
