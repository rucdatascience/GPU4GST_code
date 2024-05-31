cmake_minimum_required(VERSION 3.0.0)
project(cpu_code VERSION 0.1.0 LANGUAGES C CXX)


aux_source_directory(CPU_code DIR_SRCS)

include_directories(rucgraph CPU_code)

add_executable(test CPU_code/try.cpp)


