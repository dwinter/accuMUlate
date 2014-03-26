cmake_minimum_required(VERSION 2.6)
project(accuMUlate)

set(CMAKE_CXX_FLAGS  "-std=c++11")
find_package( Boost COMPONENTS program_options REQUIRED )

include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories(${Boost_INCLUDE_DIR})

add_executable(accuMUlate main.cc model.cc parsers.cc)
target_link_libraries(accuMUlate bamtools bamtools-utils ${Boost_LIBRARIES})

#add_executable(filter_vars post_processor.cc parsers.cc)
#target_link_libraries(filter_vars bamtools bamtools-utils ${Boost_LIBRARIES})

