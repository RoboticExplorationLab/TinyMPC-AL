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

  tiny_lqr_lti.h
  tiny_lqr_lti.c

  tiny_lqr_ltv.h
  tiny_lqr_ltv.c

  tiny_mpc_lti.h
  tiny_mpc_lti.c

  tiny_mpc_ltv.h
  tiny_mpc_ltv.c

  tiny_utils.h
  tiny_utils.c

  tiny_struct.h
  tiny_struct.c

  tiny_cost.h
  tiny_cost.c

  tiny_dynamics_lti.h
  tiny_dynamics_lti.c

  tiny_dynamics_ltv.h
  tiny_dynamics_ltv.c

  tiny_constraint.h
  tiny_constraint.c

  bicycle_5d.h
  bicycle_5d.c

  planar_quadrotor.h
  planar_quadrotor.c
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

target_include_directories(tinympc
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )

  add_library(tinympc::tinympc ALIAS tinympc)