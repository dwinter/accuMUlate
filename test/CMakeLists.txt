
################################
# GTest
################################

enable_testing()
include( CTest )
#include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

include_directories(${GTEST_INCLUDE_DIRS})
#
#    add_executable(foo foo.cc)
#    target_link_libraries(foo ${GTEST_BOTH_LIBRARIES})
#
#    add_test(AllTestsInFoo foo)

#/home/steven/Postdoc2/Project_accuMUlate/accuMUlate/include/gtest-1.7.0

################################
# Unit Tests
################################
# Add test cpp file
#SET(TEST_SOURCE_FILES ${TEST_SOURCE_FILES} testA.cc testB.cc)
SET(TEST_SOURCE_FILES ${TEST_SOURCE_FILES}

    ${CMAKE_CURRENT_SOURCE_DIR}/boost_input_utils_test.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/data_struct_test.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/model_test.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/parser_test.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/stats_test.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/unittest_utils.cc
         )


SET(LOCAL_LIBS LibDistributions LibIoData LibMutations LibStats)


add_executable( UnitTests ${TEST_SOURCE_FILES} ${SOURCE_FILES} )
# Link test executable against gtest & gtest_main
# target_link_libraries(UnitTests gtest gtest_main ${LIBS} ${LOCAL_LIBS})
target_link_libraries(UnitTests ${GTEST_BOTH_LIBRARIES} ${LIBS} ${LOCAL_LIBS} ${CMAKE_THREAD_LIBS_INIT} )
add_test( UnitTestsALL UnitTests )
##
## cd build/test
## ctest


