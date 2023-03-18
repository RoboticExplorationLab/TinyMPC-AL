function(add_example name)
  set(EXAMPLE_NAME ${name}_example)
  add_executable(${EXAMPLE_NAME}
    ${EXAMPLE_NAME}.c
    ${ARGN}
    )
  target_link_libraries(${EXAMPLE_NAME}
    PRIVATE
    slap::slap
    tinympc::tinympc
    gtest::gtest
    )
endfunction()

add_example(bicycle bicycle_5d.h bicycle_5d.c)