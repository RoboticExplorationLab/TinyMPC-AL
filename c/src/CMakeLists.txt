set(BUILD_SHARED_LIBS NO)

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
)

target_link_libraries(tinympc
  PUBLIC
  slap::slap
  )

target_include_directories(tinympc
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)

add_library(tinympc::tinympc ALIAS tinympc)