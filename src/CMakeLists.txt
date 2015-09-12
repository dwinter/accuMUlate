
if(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
  message("Using /usr/local/lib/ on ${CMAKE_SYSTEM_NAME}")
  link_directories(/usr/local/lib)
endif()


#ADD_SUBDIRECTORY(algorithm)
ADD_SUBDIRECTORY(distributions)
#ADD_SUBDIRECTORY(evolution_models)
ADD_SUBDIRECTORY(io_data)
ADD_SUBDIRECTORY(mutations)

