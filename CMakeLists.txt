cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

message("generator:${CMAKE_GENERATOR}")
message("C++ compiler:${CMAKE_CXX_COMPILER}")
project(DPBF VERSION 0.0.1 LANGUAGES CXX)
enable_language(CUDA)
find_package(CUDA REQUIRED)
# * Determine GPU architectures
if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CUDA_STANDARD 17)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")


# * Enable the CMake CUDA language

add_executable(${PROJECT_NAME} src/try_gpu.cpp src/DPQ.cu)
target_include_directories(${PROJECT_NAME} PUBLIC ${CUDA_INCLUDE_DIRS})
set(BOOST_ROOT "/home/lijiayu/boost")
# 禁用搜索系统路径
set(Boost_NO_SYSTEM_PATHS ON) # 非常重要，否则会导致找不到 Boost 库
# 查找Boost库，指定需要的组件
find_package(Boost REQUIRED)
message("${Boost_INCLUDE_DIRS}")
message("${Boost_LIBRARIES}")
# 包含Boost头文件目录
include_directories(${Boost_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC rucgraph)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)
aux_source_directory(src DIR_SRCS)
set(SPECIAL_SOURCE_FILE src/DPQ.cu)
 

# 为特定源文件设置特殊的编译命令
set_source_files_properties(
    ${SPECIAL_SOURCE_FILE}
    PROPERTIES
    COMPILE_FLAGS "--expt-extended-lambda --expt-relaxed-constexpr"
)
