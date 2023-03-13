if (SLAP_CLANG_TIDY AND (CMAKE_C_COMPILER_ID STREQUAL "Clang"))
  message(STATUS "Running clang-tidy on source directory.")
  set(CMAKE_C_CLANG_TIDY
    clang-tidy;
    -export-fixes=tidy-fixes.yaml;
    -header-filter=.;
    --config-file=${PROJECT_SOURCE_DIR}/.clang-tidy;
    )
endif ()

add_library(slap
  matrix.h
  matrix.c

  errors.h
  errors.c

  matrix_checks.h
  matrix_checks.c

  printing.h
  printing.c

  new_matrix.h
  new_matrix.c

  copy_matrix.h
  copy_matrix.c

  unary_ops.h
  unary_ops.c

  binary_ops.h
  binary_ops.c

  vector_ops.h
  vector_ops.c

  vector_products.h
  vector_products.c

  iterator.h
  iterator.c

  strided_matrix.h
  strided_matrix.c

  function_mapping.h
  function_mapping.c

  matmul.h
  matmul.c

  cholesky.h
  cholesky.c qr.c qr.h)

target_include_directories(slap
  INTERFACE
  $<BUILD_INTERFACE:${slap_SOURCE_DIR}
  $<INSTALL_INTERFACE:"include">
  )

# Link math library
target_link_libraries(slap PUBLIC m)

add_library(slap::slap ALIAS slap)