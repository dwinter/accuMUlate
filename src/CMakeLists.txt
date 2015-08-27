

SET(SOURCE_FILES CACHE INTERNAL "Global")
#message("A: ${SOURCE_FILES}")

#file(GLOB helloworld_SRC
#    "*.h"
#    "*.cpp"
#)
#add_executable(helloworld ${helloworld_SRC})

#SET(SOURCE_FILES ${SOURCE_FILES}
#    ${CMAKE_CURRENT_SOURCE_DIR}/constant.h
#    ${CMAKE_CURRENT_SOURCE_DIR}/lookup.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/setup_utils.cpp
#
#    CACHE INTERNAL "Global"
#    )

#subdirs(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#list (APPEND SOURCE_LIST ${SOURCE_FILES})
#file(GLOB_RECURSE SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

#message("BB: ${SOURCE_FILES} ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ")

if(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
  message("Using /usr/local/lib/ on ${CMAKE_SYSTEM_NAME}")
  link_directories(/usr/local/lib)
endif()


#ADD_SUBDIRECTORY(algorithm)
#ADD_SUBDIRECTORY(distributions)
#ADD_SUBDIRECTORY(evolution_models)
ADD_SUBDIRECTORY(io_data)
#ADD_SUBDIRECTORY(mutations)


#SET(LOCAL_LIBS LibEvolutionModels LibDistributions LibIoData LibMutations)

#SET(LOCAL_LIBS LibIoData )

#message("C: ${SOURCE_FILES} ")
#message(${CMAKE_CURRENT_SOURCE_DIR})

#add_executable(accuMUlate main.cc ${SOURCE_FILES})
#target_link_libraries(accuMUlate ${LIBS} ${LOCAL_LIBS})
#
#
#add_executable(MainEM main_em.cc ${SOURCE_FILES} )
#
#target_link_libraries(MainEM ${LIBS} ${LOCAL_LIBS})
##target_link_libraries(MainEM profiler)
#
#add_executable(AlignmentPreprocessor main_alignment_preprocessor.cc ${SOURCE_FILES})
#target_link_libraries(AlignmentPreprocessor ${LIBS} ${LOCAL_LIBS})
##target_link_libraries(AlignmentPreprocessor profiler)
#
#add_executable(MainRefitData main_refit_data.cc ${SOURCE_FILES} )
#target_link_libraries(MainRefitData ${LIBS} ${LOCAL_LIBS})
#
#
#
##add_executable(pp utils/post_processor.cc parsers.cc model.cc)
##target_link_libraries(pp ${Boost_LIBRARIES} ${Bamtools_LIBRARIES})
##link_directories("${Bamtools_LIBRARY_DIR}/bamtools/")
##target_link_libraries(pp ${Boost_LIBRARIES} bamtools bamtools-utils )
#
#
##add_executable(freq_seq utils/freq_seq.cc parsers.cc)
##target_link_libraries(pp bamtools bamtools-utils)
