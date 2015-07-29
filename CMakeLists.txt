cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(accuMUlate)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS  "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -O3")



########## FIND_PACKAGE ##########
#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/Modules/")
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_INSTALL_PREFIX} ~/lib ~/include ~)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakeModule/")
find_package( Boost COMPONENTS program_options thread system REQUIRED )
find_package( Eigen3 3 REQUIRED )
find_package( Bamtools REQUIRED )

FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
ENDIF(ZLIB_FOUND)

set(LIBS ${LIBS} ${Boost_LIBRARIES} ${Bamtools_LIBRARIES} )

set(INCLUDE_DIR ${INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR}  ${EIGEN3_INCLUDE_DIR})
message(${Bamtools_INCLUDE_DIR})
message(${Bamtools_LIBRARIES})

include_directories(SYSTEM ${INCLUDE_DIR})

include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories(${Boost_INCLUDE_DIR})
include_directories("./")

add_executable(accuMUlate main.cc model.cc parsers.cc)
target_link_libraries(accuMUlate ${LIBS})

#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${LIBS})

#add_executable(denom utils/denom.cc parsers.cc model.cc)
#target_link_libraries(denom ${LIBS})