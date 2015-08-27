cmake_minimum_required(VERSION 2.6)
project(accuMUlate)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from
 there. You may need to remove CMakeCache.txt. ")
endif()



SET(CMAKE_CXX_FLAGS  "-std=c++11")

SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ~/lib ~/include ${CMAKE_INSTALL_PREFIX} ~)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/Modules/")
find_package( Boost COMPONENTS program_options REQUIRED )
find_package( Bamtools REQUIRED )
find_package( Eigen3 3 REQUIRED )


set(LIBS ${LIBS} ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
#set(LIBS ${LIBS} ${Boost_LIBRARIES},  "${CMAKE_SOURCE_DIR}/third-party/bamtools/lib") try to delete this completely
set(INCLUDE_DIR ${INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR}  ${EIGEN3_INCLUDE_DIR})

include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories(SYSTEM ${INCLUDE_DIR})
include_directories("./")

#include_directories(${CMAKE_SOURCE_DIR}/src) not at moment, maybe add it later





add_executable(accuMUlate main.cc model.cc boost_input_utils.cc parsers.cc)
target_link_libraries(accuMUlate ${LIBS})

#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${LIBS})

#add_executable(denom utils/denom.cc parsers.cc model.cc)
#target_link_libraries(denom ${LIBS})
