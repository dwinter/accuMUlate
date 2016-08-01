cmake_minimum_required(VERSION 2.6)
project(accuMUlate)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from
 there. You may need to remove CMakeCache.txt. ")
endif()



SET(CMAKE_CXX_FLAGS  "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE  "-std=c++11 -O3")
set(CMAKE_CXX_FLAGS_RelWithDebInfo "-std=c++11 -O3 -Wall -Wextra -pedantic -Wno-deprecated -Wno-deprecated-declarations -Wno-deprecated")
set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -g  -O0 -Wvla -Wall  -Wextra -pedantic ")



SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ~/local/lib ~/local/include ${CMAKE_INSTALL_PREFIX} ~)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/Modules/")
find_package( Boost COMPONENTS program_options REQUIRED )
find_package( Bamtools REQUIRED )
#find_package( Eigen3 3 REQUIRED )


set(LIBS ${LIBS} ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
#set(LIBS ${LIBS} ${Boost_LIBRARIES},  "${CMAKE_SOURCE_DIR}/third-party/bamtools/lib") try to delete this completely
set(INCLUDE_DIR ${INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR})#  ${EIGEN3_INCLUDE_DIR})

#include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
#include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories(SYSTEM ${INCLUDE_DIR})
include_directories("./")

#include_directories(${CMAKE_SOURCE_DIR}/src) not at moment, maybe add it later

SET(SOURCE_FILES CACHE INTERNAL "Global")
SET(SOURCE_FILES ${SOURCE_FILES}
        ${CMAKE_CURRENT_SOURCE_DIR}/boost_input_utils.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/data_struct.h
        ${CMAKE_CURRENT_SOURCE_DIR}/model.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/parsers.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/variant_visitor.cc
        CACHE INTERNAL "Global"
        )

ADD_SUBDIRECTORY(src)
SET(LOCAL_LIBS LibDistributions LibIoData LibMutations LibStats)
#SET(LOCAL_LIBS LibEvolutionModels LibDistributions LibIoData LibMutations)



add_executable(accuMUlate main.cc ${SOURCE_FILES})
target_link_libraries(accuMUlate ${LIBS} ${LOCAL_LIBS})

#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${LIBS})

#add_executable(denom utils/denom.cc parsers.cc model.cc)
#target_link_libraries(denom ${LIBS})



find_package (GTest)
IF(GTEST_FOUND)

    #1
    #ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/include/googletest-trunk/")
    #SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/include/googletest-trunk")

    find_package(Threads)

    IF(THREADS_FOUND)
        SET(GTEST_THREAD_LIB "-lpthread")
#        message("${ZLIB_INCLUDE_DIR}  ${ZLIB_LIBRARIES}  ${ZLIB_FOUND}")
#        message("${GTEST_THREAD_LIB}A${CMAKE_THREAD_LIBS_INIT}B${Threads}====")
    ENDIF(THREADS_FOUND)


    message("GTEST FOUND! Build UnitTest")
    ADD_SUBDIRECTORY(test)


    #2 Some link error??
    #SET(GTEST_ROOT "${CMAKE_SOURCE_DIR}/include/googletest-trunk")
    #find_package( GTest  )
    #IF(GTEST_FOUND)
    #    message("Found GTEST at ${GTEST_INCLUDE_DIRS} ${gtest_SOURCE_DIR} ")
    #    include_directories(${GTEST_INCLUDE_DIRS})
    #    ADD_SUBDIRECTORY(test)
    #ENDIF(GTEST_FOUND)

    #3 SVN checkout version
    #ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/gtest")

    # Defines the following variables:
    #
    #   GTEST_FOUND - Found the Google Testing framework
    #   GTEST_INCLUDE_DIRS - Include directories
    #
    # Also defines the library variables below as normal
    # variables.  These contain debug/optimized keywords when
    # a debugging library is found.
    #
    #   GTEST_BOTH_LIBRARIES - Both libgtest & libgtest-main
    #   GTEST_LIBRARIES - libgtest
    #   GTEST_MAIN_LIBRARIES - libgtest-main
    #
ENDIF(GTEST_FOUND)
