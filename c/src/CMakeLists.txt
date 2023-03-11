add_library(tinympc
  # constrained_ilqr.h
  # constrained_ilqr.c

  constrained_lqr.h
  constrained_lqr.c

  unconstrained_lqr.h
  unconstrained_lqr.c

  tiny_utils.h

  tiny_struct.h
  tiny_struct.c

  bicycle.h
  bicycle.c

  planar_quadrotor.h
  planar_quadrotor.c
  )

target_link_libraries(tinympc
  PRIVATE
  slap
  )

target_include_directories(tinympc
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )