# Made by Sam Schoedel on 3/12/2023
# I don't know the best way to include external projects using CMake yet
# Need to ask someone knowledgeable (probably Brian)
# For now just creating a new library for each external project and then 
# doing target_link_libraries(executable_name PRIVATE lib_name)

# I also don't know if this is the correct way to allow include headers to
# not require this folder name (currently 'lib') in front of the library path
# for ex #include "cjson/cJSON.h" instead of #include "lib/cjson/cJSON.h"
# include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(cjson)
add_subdirectory(slap)
