set ( PROBES_TO_GENES_SRC MapProbesToGenes.cpp affymetrix_utils.cpp )
set ( SAILFISH_MAIN_SRCS 
#affymetrix_utils.cpp 
#nnls.cpp
#G12.cpp
#myMath.cpp 
#H12.cpp
#Indexer.cpp
PerfectHashIndexer.cpp
IndexedCounter.cpp
sailfish.cpp 
)

set( ANALYZE_TRANSCRIPT_GRAPH 
AnalyzeTranscriptGraph.cpp
IndexedCounter.cpp
PerfectHashIndexer.cpp
)

set( BUILD_LUT_SRCS 
BuildLUT.cpp
IndexedCounter.cpp
PerfectHashIndexer.cpp
)

set( COMPUTE_TGRAPH_SRCS
ComputeTranscriptGraph.cpp
IndexedCounter.cpp
PerfectHashIndexer.cpp
)


set ( BUILD_INDEX_SRCS Indexer.cpp )
set ( TEST_SRCS tsnnls_is_stupid.cpp )

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/nnls
${GAT_SOURCE_DIR}/external/tsnnls
${GAT_SOURCE_DIR}/external/mlrmcl
${GAT_SOURCE_DIR}/external/g2log/src
${ZLIB_INCLUDE_DIR}
${SEQAN_PATH}/core/include
${SEQAN_PATH}/extra/include
/home/robp/dist/include
/home/robp/dist/include/jellyfish-1.1.6
/opt/local/include/jellyfish-1.1.6
/fs/szattic-asmg6/rob/dist/include/jellyfish-1.1.6
/fs/szattic-asmg6/rob/dist/include
/fs/szattic-asmg6/rob/dist/include
/opt/local/include
/opt/local/include/libxml2
/opt/local/include/eigen3
)

link_directories( 
${GAT_SOURCE_DIR}/lib
/home/robp/dist/lib
/home/robp/dist/lib64
/opt/local/lib
/opt/local/lib/gcc47
/fs/szattic-asmg6/rob/dist/lib64
/fs/szattic-asmg6/rob/dist/lib
)

if (Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
endif()

if (ZLIB_FOUND)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_HAS_ZLIB")
endif (ZLIB_FOUND)

#add_executable( MapProbesToGenes ${PROBES_TO_GENES_SRC} )
#target_link_libraries( MapProbesToGenes ${Boost_LIBRARIES} ${ZLIB_LIBRARY} )

message ( "Boost Libraries are ${Boost_LIBRARIES}" )
message ( "Sailfish sources are ${SAILFISH_MAIN_SRCS}")
add_executable( Sailfish ${SAILFISH_MAIN_SRCS} )
target_link_libraries( Sailfish 
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish pthread 
	#cppa mlrmcl 
	nnls
	#tsnnls lapack
    gomp lib_activeobject lib_g2logger blas m tbb)#tsnnls shotgun )

add_executable( BuildLUT ${BUILD_LUT_SRCS} )
target_link_libraries( BuildLUT 
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish pthread 
    gomp lib_activeobject lib_g2logger blas m tbb)

add_executable( ComputeTranscriptGraph ${COMPUTE_TGRAPH_SRCS} )
target_link_libraries( ComputeTranscriptGraph
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish pthread 
    gomp lib_activeobject lib_g2logger blas m tbb nnls)


add_executable( AnalyzeTranscriptGraph ${ANALYZE_TRANSCRIPT_GRAPH} )
target_link_libraries( AnalyzeTranscriptGraph 
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish pthread 
    gomp lib_activeobject lib_g2logger blas m tbb)


#add_executable( BuildIndex ${BUILD_INDEX_SRCS} )
#target_link_libraries( BuildIndex ${Boost_LIBRARIES} jellyfish pthread m tbb)


#add_executable( TestTsnnls ${TEST_SRCS} )
#target_link_libraries( TestTsnnls ${Boost_LIBRARIES} ${ZLIB_LIBRARY} blas lapack gomp lib_activeobject lib_g2logger tsnnls )

