add_library(LibIoData
    local_bamtools/bamtools_fasta.cpp
    local_bamtools/bamtools_pileup_engine.cpp

)
#target_link_libraries(LibIoData ${LIBS} )
target_link_libraries(LibIoData ${Bamtools_LIBRARIES})

