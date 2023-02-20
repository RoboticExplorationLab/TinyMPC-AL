if (SLAP_CLANG_TIDY AND (CMAKE_C_COMPILER_ID STREQUAL "Clang")) 
  message(STATUS "Running clang-tidy on source directory.")
  set(CMAKE_C_CLANG_TIDY
    clang-tidy;
    -export-fixes=tidy-fixes.yaml;
    -header-filter=.;
    --config-file=${PROJECT_SOURCE_DIR}/.clang-tidy;
  )
endif()

add_library(slap
    matrix.h
    matrix.c

    linalg.h
    linalg.c

    submatrix.h
    submatrix.c
)

## Places header files in a "slap" subfolder in the build tree
##   and adds it to the include directories
target_include_directories(
  slap INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:"include">
)
#target_include_directories(slap INTERFACE
#    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
#    $<INSTALL_INTERFACE:"include">
#)
#file(COPY ${CMAKE_CURRENT_LIST_DIR}/
#  DESTINATION ${PROJECT_BINARY_DIR}/slap
#  FILES_MATCHING
#    PATTERN *.h
#)

# add_target_to_install(slap)
