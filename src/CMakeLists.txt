set ( SAILFISH_MAIN_SRCS 
Sailfish.cpp
)

# Build the specific Sailfish commands as a 
# shared libraray (for some reason, building
# them into a single executable was causing 
# multiple symbold definition errors).
set (SAILFISH_LIB_SRCS
QuantificationDriver.cpp
PerfectHashIndexer.cpp
BuildLUT.cpp
IndexedCounter.cpp
GenomicFeature.cpp
JellyfishMerCounter.cpp
VersionChecker.cpp
Utils.cpp
)

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/jellyfish-1.1.10
${ZLIB_INCLUDE_DIR}
${Boost_INCLUDE_DIRS}
)

link_directories( 
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
${Boost_LIBRARY_DIRS}
)

# Build the Sailfish library
add_library(sailfish SHARED ${SAILFISH_LIB_SRCS})
# Build the Sailfish executable
add_executable( Sailfish ${SAILFISH_MAIN_SRCS} GenomicFeature.cpp )
# Link the executable
target_link_libraries( Sailfish 
	sailfish 
	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
	cmph # perfect hashing library
	jellyfish-1.1 pthread 
    gomp m tbb
)


set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

install(DIRECTORY 
        ${GAT_SOURCE_DIR}/external/install/lib/ 
        DESTINATION ${INSTALL_LIB_DIR}
    )

# install(FILES ${Boost_LIBRARIES}
# 	           DESTINATION ${INSTALL_LIB_DIR})

install(TARGETS Sailfish sailfish
                RUNTIME DESTINATION bin 
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
        )

# install(FILES ${GAT_SOURCE_DIR}/scripts/SFPipeline.py DESTINATION scripts/SFPipeline.py )
# install(FILES ${GAT_SOURCE_DIR}/experimental_configs/SRX016368_25mers.cfg DESTINATION experimental_configs/SRX016368_25mers.cfg )

include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "Sailfish")
set(CPACK_SOURCE_PACKAGE_NAME "Sailfish-Src")
set(CPACK_PACKAGE_VENDOR "Carnegie Mellon University")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Sailfish - Ultrafast alignment-free RNA-seq quantification")
set(CPACK_PACKAGE_VERSION "0.1.0")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")

####
#
# Deprecated or currently unused
#
####


# # use, i.e. don't skip the full RPATH for the build tree
# SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# # when building, don't use the install RPATH already
# # (but later on when installing)
# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
# MESSAGE( "LINK DIRECTORIES : ${CMAKE_LIBRARY_PATH}")
# SET(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_PATH}")

# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# # the RPATH to be used when installing, but only if it's not a system directory
# LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# IF("${isSystemDir}" STREQUAL "-1")
#    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# ENDIF("${isSystemDir}" STREQUAL "-1")


# set( BUILD_LUT_SRCS 
# BuildLUT.cpp
# IndexedCounter.cpp
# PerfectHashIndexer.cpp
# )

# add_executable( BuildLUT ${BUILD_LUT_SRCS} )
# target_link_libraries( BuildLUT 
# 	${Boost_LIBRARIES} ${ZLIB_LIBRARY} 
# 	cmph # perfect hashing library
# 	jellyfish-1.1 pthread 
#     gomp lib_activeobject lib_g2logger m tbb)

# INSTALL(TARGETS BuildLUT DESTINATION ${PROJECT_SOURCE_DIR}/bin COMPONENT comp_buildlut)

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

