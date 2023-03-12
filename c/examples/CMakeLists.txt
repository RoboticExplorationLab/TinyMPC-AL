function(add_example name)
  set(TEST_NAME ${name}_example)
  add_executable(${TEST_NAME}
    ${TEST_NAME}.c
    )
  target_link_libraries(${TEST_NAME}
    PRIVATE
    simpletest
    slap
    tinympc
    )
endfunction()

add_example(bicycle)