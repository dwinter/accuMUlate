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
find_package( Eigen3 3 REQUIRED )


set(LIBS ${LIBS} ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
set(INCLUDE_DIR ${INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR}  ${EIGEN3_INCLUDE_DIR})

include_directories(SYSTEM ${INCLUDE_DIR})
include_directories("./")


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



add_executable(accuMUlate main.cc ${SOURCE_FILES})
target_link_libraries(accuMUlate ${LIBS} ${LOCAL_LIBS})

add_executable(denominate denom.cc ${SOURCE_FILES})
target_link_libraries(denominate ${LIBS} ${LOCAL_LIBS})



find_package (GTest)
IF(GTEST_FOUND)
    
    find_package(Threads)

    IF(THREADS_FOUND)
        SET(GTEST_THREAD_LIB "-lpthread")
    ENDIF(THREADS_FOUND)


    message("GTEST FOUND! Build UnitTest")
    ADD_SUBDIRECTORY(test)

ENDIF(GTEST_FOUND)
