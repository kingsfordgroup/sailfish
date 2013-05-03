set ( SAILFISH_MAIN_SRCS 
PerfectHashIndexer.cpp
IndexedCounter.cpp
Sailfish.cpp 
)

set( BUILD_LUT_SRCS 
BuildLUT.cpp
IndexedCounter.cpp
PerfectHashIndexer.cpp
)

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-1.1.10
${GAT_SOURCE_DIR}/external/g2log/src
${ZLIB_INCLUDE_DIR}
${Boost_INCLUDE_DIRS}
)

link_directories( 
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
${Boost_LIBRARY_DIRS}
)

add_executable( Sailfish ${SAILFISH_MAIN_SRCS} )
target_link_libraries( Sailfish 
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish-1.1 pthread 
    gomp lib_activeobject lib_g2logger m tbb)

add_executable( BuildLUT ${BUILD_LUT_SRCS} )
target_link_libraries( BuildLUT 
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish-1.1 pthread 
    gomp lib_activeobject lib_g2logger m tbb)

INSTALL(TARGETS Sailfish DESTINATION ${PROJECT_SOURCE_DIR}/bin COMPONENT comp_sailfish)
INSTALL(TARGETS BuildLUT DESTINATION ${PROJECT_SOURCE_DIR}/bin COMPONENT comp_buildlut)

# set ( BUILD_INDEX_SRCS Indexer.cpp )

# set( ANALYZE_TRANSCRIPT_GRAPH 
# AnalyzeTranscriptGraph.cpp
# IndexedCounter.cpp
# PerfectHashIndexer.cpp
# )

# set ( PROBES_TO_GENES_SRC MapProbesToGenes.cpp affymetrix_utils.cpp )
# set( COMPUTE_TGRAPH_SRCS
# ComputeTranscriptGraph.cpp
# IndexedCounter.cpp
# PerfectHashIndexer.cpp
# )

# add_executable( ComputeTranscriptGraph ${COMPUTE_TGRAPH_SRCS} )
# target_link_libraries( ComputeTranscriptGraph
# 	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
# 	cmph # perfect hashing library
# 	jellyfish pthread 
#     gomp lib_activeobject lib_g2logger blas m tbb nnls)


# add_executable( AnalyzeTranscriptGraph ${ANALYZE_TRANSCRIPT_GRAPH} )
# target_link_libraries( AnalyzeTranscriptGraph 
# 	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
# 	cmph # perfect hashing library
# 	jellyfish pthread 
#     gomp lib_activeobject lib_g2logger blas m tbb)


#add_executable( BuildIndex ${BUILD_INDEX_SRCS} )
#target_link_libraries( BuildIndex ${Boost_LIBRARIES} jellyfish pthread m tbb)


#add_executable( TestTsnnls ${TEST_SRCS} )
#target_link_libraries( TestTsnnls ${Boost_LIBRARIES} ${ZLIB_LIBRARY} blas lapack gomp lib_activeobject lib_g2logger tsnnls )

