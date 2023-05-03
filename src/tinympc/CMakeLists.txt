# if (CLANG_TIDY AND (CMAKE_C_COMPILER_ID STREQUAL "Clang"))
#   message(STATUS "Running clang-tidy on source directory.")
#   set(CMAKE_C_CLANG_TIDY
#     clang-tidy;
#     -export-fixes=tidy-fixes.yaml;
#     -header-filter=.;
#     --config-file=${PROJECT_SOURCE_DIR}/.clang-tidy;
#     )
# endif ()

add_library(tinympc
  constants.h 
  
  utils.h
  utils.c

  types.h

  model.h
  model.c

  auxil.h
  auxil.c
  
  cost_lqr.h
  cost_lqr.c

  lqr.h
  lqr.c

  # lqr_ltv.h
  # lqr_ltv.c

  constraint_linear.h
  constraint_linear.c

  al_lqr.h
  al_lqr.c

  # mpc_lti.h
  # mpc_lti.c

  # mpc_ltv.h
  # mpc_ltv.c



  # dynamics_lti.h
  # dynamics_lti.c

  # dynamics_ltv.h
  # dynamics_ltv.c


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