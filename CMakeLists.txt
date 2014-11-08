cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(accuMUlate)

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()



set(CMAKE_CXX_FLAGS  "-std=c++11")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakeModule/")

find_package( Boost COMPONENTS program_options REQUIRED )
find_package( Eigen3 3 REQUIRED )
find_package( Bamtools REQUIRED )
find_package( gsl REQUIRED )

FIND_PACKAGE(ZLIB)
IF(ZLIB_FOUND)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
ENDIF(ZLIB_FOUND)


set(LIBS ${LIBS} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${Bamtools_LIBRARIES})
set(INCLUDE_DIR ${INCLUDE_DIR} ${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIR} ${Bamtools_INCLUDE_DIR}  ${EIGEN3_INCLUDE_DIR})


include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories(${INCLUDE_DIR})

ADD_SUBDIRECTORY(distributions)


add_executable(accuMUlate main.cc model.cc parsers.cc)
target_link_libraries(accuMUlate ${LIBS} )
#target_link_libraries(accuMUlate distributions)
#target_link_libraries(accuMUlate bamtools bamtools-utils ${LIBS} distributions)


#add_executable(pp utils/post_processor.cc parsers.cc model.cc)
#target_link_libraries(pp ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
#link_directories("${Bamtools_LIBRARY_DIR}/bamtools/")
#target_link_libraries(pp ${Boost_LIBRARIES} bamtools bamtools-utils )


#add_executable(freq_seq utils/freq_seq.cc parsers.cc)
#target_link_libraries(pp bamtools bamtools-utils)


IF(USER_DEBUG)
	MESSAGE("===== inc ${INCLUDE_DIR}
	${Boost_INCLUDE_DIR} 
	${BamTools_INCLUDE_DIR}
	")
	
	message("===== libs ${LIBS}
	${Boost_LIBRARIES} 
	${GSL_LIBRARIES} 
	${Bamtools_LIBRARIES}
	${Bamtools_LIBRARY_DIR}
	")
ENDIF (USER_DEBUG)	