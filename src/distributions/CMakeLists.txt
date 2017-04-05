#add_executable(accuMUlate DirichletMultinomialDistribution.cpp)


add_library(LibDistributions DirichletMultinomialDistribution.cpp)


SET(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/DirichletMultinomialDistribution.cpp
    #PARENT_SCOPE
    CACHE INTERNAL "Global"
    )
