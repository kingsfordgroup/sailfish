set ( PROBES_TO_GENES_SRC MapProbesToGenes.cpp affymetrix_utils.cpp )
set ( SAILFISH_MAIN_SRCS sailfish.cpp affymetrix_utils.cpp )
set ( TEST_SRCS tsnnls_is_stupid.cpp )

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/mlrmcl
${GAT_SOURCE_DIR}/external/g2log/src
${ZLIB_INCLUDE_DIR}
${SEQAN_PATH}/core/include
${SEQAN_PATH}/extra/include
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
/fs/szattic-asmg6/rob/dist/lib64
/fs/szattic-asmg6/rob/dist/lib
/opt/local/lib
/opt/local/lib/gcc47
)

if (Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
endif()

if (ZLIB_FOUND)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_HAS_ZLIB")
endif (ZLIB_FOUND)

add_executable( MapProbesToGenes ${PROBES_TO_GENES_SRC} )
message ( "Boost Libraries are ${Boost_LIBRARIES}" )
target_link_libraries( MapProbesToGenes ${Boost_LIBRARIES} ${ZLIB_LIBRARY} )

add_executable( Sailfish ${SAILFISH_MAIN_SRCS} )
target_link_libraries( Sailfish ${Boost_LIBRARIES} ${ZLIB_LIBRARY} jellyfish pthread cppa mlrmcl gomp lib_activeobject lib_g2logger 
nnls blas lapack m )#tsnnls shotgun )

#add_executable( TestTsnnls ${TEST_SRCS} )
#target_link_libraries( TestTsnnls ${Boost_LIBRARIES} ${ZLIB_LIBRARY} blas lapack gomp lib_activeobject lib_g2logger tsnnls )

