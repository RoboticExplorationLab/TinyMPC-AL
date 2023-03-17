# TODO (sschoedel): add support for GNU build
if (TINY_CLANG_TIDY AND (CMAKE_C_COMPILER_ID STREQUAL "Clang"))
  message(STATUS "Running clang-tidy on source directory.")
  set(CMAKE_C_CLANG_TIDY
    clang-tidy;
    -export-fixes=tidy-fixes.yaml;
    -header-filter=.;
    --config-file=${PROJECT_SOURCE_DIR}/.clang-tidy;
    )
endif ()

add_library(tinympc
  # constrained_ilqr.h
  # constrained_ilqr.c

  lqr_lti.h
  lqr_lti.c

  lqr_ltv.h
  lqr_ltv.c

  mpc_lti.h
  mpc_lti.c

  mpc_ltv.h
  mpc_ltv.c

  utils.h
  utils.c

  data_struct.h
  data_struct.c

  cost_lqr.h
  cost_lqr.c

  dynamics_lti.h
  dynamics_lti.c

  dynamics_ltv.h
  dynamics_ltv.c

  constraint_linear.h
  constraint_linear.c
  )

target_include_directories(tinympc
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:"include">
  )

target_link_libraries(tinympc
  PRIVATE
  slap::slap
  )

add_library(tinympc::tinympc ALIAS tinympc)