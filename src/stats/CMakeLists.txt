SET(SOURCE_FILES ${SOURCE_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/stats.cc
    #PARENT_SCOPE
    CACHE INTERNAL "Global"
)


add_library(LibStats
    stats.cc
)

