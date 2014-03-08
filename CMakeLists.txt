cmake_minimum_required(VERSION 2.6)
project(accuMUlate)

set(CMAKE_CXX_FLAGS  "-std=c++11")
include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src"
    "${CMAKE_SOURCE_DIR}/third-party/")
add_executable(run_tetma parse_bam.cc model.cc)
target_link_libraries(run_tetma bamtools bamtools-utils)
#target_link_librari
