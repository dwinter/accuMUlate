cmake_minimum_required(VERSION 2.6)
project(accuMUlate)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/Modules/")
set(CMAKE_CXX_FLAGS  "-std=c++11")
find_package( Boost COMPONENTS program_options REQUIRED )
find_package( Bamtools REQUIRED )



set(LIBS ${LIBS} ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
#set(LIBS ${LIBS} ${Boost_LIBRARIES},  "${CMAKE_SOURCE_DIR}/third-party/bamtools/lib")
include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories(${Boost_INCLUDE_DIR})
include_directories("./")

add_executable(accuMUlate main.cc model.cc boost_input_utils.cc parsers.cc)
target_link_libraries(accuMUlate ${LIBS})

#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${LIBS})

#add_executable(denom utils/denom.cc parsers.cc model.cc)
#target_link_libraries(denom ${LIBS})
