add_library(tinympc
  riccati.h
  riccati.c

  util.h

  tiny_struct.h
  tiny_struct.c

  augmented_lagrangian_lqr.h
  augmented_lagrangian_lqr.c
  )

target_link_libraries(tinympc
  PRIVATE
  slap
  )

target_include_directories(tinympc
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )