cmake_minimum_required(VERSION 2.6)
project(accuMUlate)

set(CMAKE_CXX_FLAGS  "-std=c++11")
find_package( Boost COMPONENTS program_options REQUIRED )

include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories(${Boost_INCLUDE_DIR})
include_directories("./")

add_executable(accuMUlate main.cc model.cc parsers.cc)
target_link_libraries(accuMUlate bamtools bamtools-utils ${Boost_LIBRARIES})

#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp bamtools bamtools-utils ${Boost_LIBRARIES})

#add_executable(freq_seq utils/freq_seq.cc parsers.cc)
#target_link_libraries(pp bamtools bamtools-utils)
