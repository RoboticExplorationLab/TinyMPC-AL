add_library(tinympc
  riccati.h
  riccati.c

  util.h

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