set ( PROBES_TO_GENES_SRC MapProbesToGenes.cpp affymetrix_utils.cpp )
set ( SAILFISH_MAIN_SRCS sailfish.cpp affymetrix_utils.cpp )

include_directories( 
${GAT_SOURCE_DIR}/include
${ZLIB_INCLUDE_DIR}
${SEQAN_PATH}/core/include
${SEQAN_PATH}/extra/include
/opt/local/include/jellyfish-1.1.6
/opt/local/include
/opt/local/include/libxml2
/opt/local/include/eigen3
)

link_directories( 
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
target_link_libraries( MapProbesToGenes ${Boost_LIBRARIES} ${ZLIB_LIBRARY})

add_executable( Sailfish ${SAILFISH_MAIN_SRCS} )
target_link_libraries( Sailfish ${Boost_LIBRARIES} ${ZLIB_LIBRARY} jellyfish )
