
aux_source_directory(. FILTER_SRCS)

MESSAGE(STATUS "building filter")

list(APPEND CMAKE_CXX_FLAGS "-std=c++11")

# set(OpenCV_DIR "/data/xygpu/software/opencv-4.6.0/release/include/opencv4")
find_package(OpenCV REQUIRED)
find_package(CUDA QUIET REQUIRED)
include_directories(${PROJECT_SOURCE_DIR} ${OpenCV_INCLUDE_DIRS})

# set(OpenCV_INCLUDE_PATH "/data/xygpu/software/opencv-4.6.0/release/include/opencv4")
# include_directories(${PROJECT_SOURCE_DIR} ${OpenCV_INCLUDE_PATH})
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-G -g -O0 -lcufftw -lcufft -lculibos)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-G -g -O0 -lcufftw -lcufft -lculibos)
list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_86,code=sm_86)
cuda_add_library (filter ${FILTER_SRCS})


# target_link_libraries(filter)
target_link_libraries(filter ${OpenCV_LIBS} ${CUDA_LIBRARIES} -L/usr/local/cuda/lib64 -lcufftw) #需要系统能够找到opencv才可打开

set_target_properties(filter PROPERTIES OUTPUT_NAME "filter")
set_target_properties(filter PROPERTIES CLEAN_DIRECT_OUTPUT 1)