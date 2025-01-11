
MESSAGE(STATUS "Building pnp_lib as a shared library")

FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h")

if(ENABLE_OMP)
       find_package(OpenMP REQUIRED)
       set(CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
       set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_C_FLAGS}")
endif(ENABLE_OMP)

find_package(MPI REQUIRED)
#find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)


set(MPI_INCLUDE_PATH "/usr/local/openmpi/include/")
include_directories(${PROJECT_SOURCE_DIR} ${MPI_INCLUDE_PATH} ${OpenCV_INCLUDE_DIRS})
# 创建一个共享库而不是可执行文件
#add_executable(pnp ${SOURCES})
add_library(pnp_lib SHARED ${SOURCES})

# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -fopenmp")
# 找到ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#target_link_libraries(pnp ${MPI_LIBRARIES} mrcmx ${OpenCV_LIBS} ${ITK_LIBRARIES})
target_link_libraries(pnp_lib ${MPI_LIBRARIES} mrcmx ${OpenCV_LIBS} ${ITK_LIBRARIES})

# add
SET(CMAKE_BUILD_TYPE "Debug") 
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

