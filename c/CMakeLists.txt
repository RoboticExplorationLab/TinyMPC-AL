cmake_minimum_required(VERSION 3.2)
list(APPEND CMAKE_MESSAGE_CONTEXT tinympc)
project(TinyMPC VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Include useful modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CMakePrintHelpers)
include(GNUInstallDirs)
# include(Functions)
include(FetchContent)

# Add CPM Dependency Manager
include(FindCPM)

# Handle default build type
set(TINY_DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Setting CMAKE_BUILD_TYPE to ${TINY_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE ${TINY_DEFAULT_BUILD_TYPE} CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

##############################
# Build
##############################

# Compile options
# TODO: Future options to optimize library for specific use of constraints
option(SET_GOAL_CONSTRAINT "Set goal constraint" OFF)
if (SET_GOAL_CONSTRAINT)
  add_definitions(-DGOAL_CONSTRAINT)
endif()
option(SET_INPUT_CONSTRAINT "Set input constraints" ON)
if (SET_INPUT_CONSTRAINT)
  add_definitions(-DINPUT_CONSTRAINT)
endif()
option(SET_STATE_CONSTRAINT "Set state constraints" ON)
if (SET_STATE_CONSTRAINT)
  add_definitions(-DSTATE_CONSTRAINT)
endif()

# Floating point precision
set(SLAP_FLOAT double CACHE STRING "Floating point precision for slap (float,double).")

add_compile_options(-Wall -Wextra -pedantic -Wno-error=unknown-pragmas)
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
add_subdirectory(examples)
add_subdirectory(src)
add_subdirectory(test)