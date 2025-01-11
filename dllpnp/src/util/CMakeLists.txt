
#set(UTIL_SRCS exception.cpp dataf.cpp)
aux_source_directory(. UTIL_SRCS)

MESSAGE(STATUS "building util")
#-fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include_directories(${PROJECT_SOURCE_DIR})

add_library (util STATIC ${UTIL_SRCS})

target_link_libraries(util)

set_target_properties(util PROPERTIES OUTPUT_NAME "util")
set_target_properties(util PROPERTIES CLEAN_DIRECT_OUTPUT 1)
