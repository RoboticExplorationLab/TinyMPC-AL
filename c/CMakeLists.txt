cmake_minimum_required(VERSION 3.0.0)
project(main)

include(FetchContent)
FetchContent_Declare(slap
    GIT_REPOSITORY https://github.com/bjack205/slap
    GIT_TAG main
)
FetchContent_MakeAvailable(slap)

add_executable(main main.c)
target_link_libraries(main PRIVATE slap::slap)