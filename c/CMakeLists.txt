cmake_minimum_required(VERSION 3.0.0)
project(tinyMPC VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Handle default build type
set(SLAP_DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Setting CMAKE_BUILD_TYPE to ${SLAP_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE ${SLAP_DEFAULT_BUILD_TYPE} CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

##############################
# Build
##############################

# Compile options
add_compile_options(-Wall -Wextra -pedantic -Wno-error=unknown-pragmas)
add_compile_options(-O2 -mavx2 -mfma)
# add_compile_options(-Wall -Wextra -Wpedantic -Werror
#                     -Wformat=2 -Wno-unused-parameter -Wshadow 
#                     -Wwrite-strings -Wstrict-prototypes -Wold-style-definition 
#                     -Wredundant-decls -Wnested-externs -Wmissing-include-dirs)

if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  add_compile_options(-Wjump-misses-init -Wlogical-op)
endif ()                    
# Make all includes relative to src/ directory
include_directories(${PROJECT_SOURCE_DIR})

# Build source files
add_subdirectory(slap)

##############################
# Examples
##############################
add_executable(getting_started examples/getting_started.c)
target_link_libraries(getting_started PRIVATE slap::slap)

add_subdirectory(examples/riccati)
