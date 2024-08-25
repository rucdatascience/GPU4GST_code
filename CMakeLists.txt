cmake_minimum_required(VERSION 3.23.1 FATAL_ERROR)
message("generator:${CMAKE_GENERATOR}")
message("C++ compiler:${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.0/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
endif()
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)
    file(DOWNLOAD https://raw.githubusercontent.com/rapidsai/rapids-cmake/branch-24.08/RAPIDS.cmake
         ${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)
endif()
include(${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)

include(rapids-cmake)
include(rapids-cpm)
include(rapids-cuda)
include(rapids-export)
include(rapids-find)
# * Determine GPU architectures
# * Enable the CMake CUDA language
project(DPBF VERSION 0.0.1 LANGUAGES CXX CUDA)
rapids_cuda_init_architectures(DPBF)
rapids_find_package(
    CUDAToolkit REQUIRED
 
)
add_executable(${PROJECT_NAME} src/try_gpu.cpp src/DPQ.cu)
target_link_libraries(${PROJECT_NAME}  CUDA::toolkit)
include_directories(rucgraph )
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)
aux_source_directory(src DIR_SRCS)
set(SPECIAL_SOURCE_FILE src/DPQ.cu)
 

# 为特定源文件设置特殊的编译命令
set_source_files_properties(
    ${SPECIAL_SOURCE_FILE}
    PROPERTIES
    COMPILE_FLAGS "--expt-extended-lambda --expt-relaxed-constexpr"
)