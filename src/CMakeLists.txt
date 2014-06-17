set ( SAILFISH_MAIN_SRCS 
Sailfish.cpp
)

#add_custom_command ( OUTPUT QSufSort.c COMMAND touch QSufSort.c)
#add_custom_command ( OUTPUT is.c COMMAND touch is.c)
#add_custom_command ( OUTPUT bwt_gen.c COMMAND touch bwt_gen.c)
#add_custom_command ( OUTPUT bwtindex.c COMMAND touch bwtindex.c)

set ( SALMON_MAIN_SRCS
QSufSort.c
is.c
bwt_gen.c
bwtindex.c
FragmentList.cpp
Salmon.cpp
BuildSalmonIndex.cpp
SalmonQuantify.cpp
)

set (BUILD_TRANSCRIPT_MAP_SRCS
BuildTranscriptMap.cpp
)

set (SALMON_READ_SRCS
FASTAParser.cpp 
SailfishStringUtils.cpp 
LibraryFormat.cpp 
ErrorModel.cpp 
FragmentLengthDistribution.cpp 
TestBAMParser.cpp 
)

# Build the specific Sailfish commands as a 
# shared libraray (for some reason, building
# them into a single executable was causing 
# multiple symbold definition errors).
set (SAILFISH_LIB_SRCS
LibraryFormat.cpp
QuantificationDriver.cpp
PerfectHashIndexer.cpp
BuildLUT.cpp
IndexedCounter.cpp
LookUpTableUtils.cpp
GenomicFeature.cpp
JellyfishMerCounter.cpp
VersionChecker.cpp
SailfishUtils.cpp
SailfishStringUtils.cpp
ComputeBiasFeatures.cpp
PerformBiasCorrection.cpp
PartitionRefiner.cpp
StreamingSequenceParser.cpp
cokus.cpp
merge_files.cc
)

include_directories( 
${GAT_SOURCE_DIR}/include
${GAT_SOURCE_DIR}/external
${GAT_SOURCE_DIR}/external/cereal/include
${GAT_SOURCE_DIR}/external/install/include
${GAT_SOURCE_DIR}/external/install/include/samtools
${GAT_SOURCE_DIR}/external/install/include/jellyfish-2.1.3
${GAT_SOURCE_DIR}/external/install/include/bwa
${GAT_SOURCE_DIR}/external/install/include/clasp
${ZLIB_INCLUDE_DIR}
${TBB_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS}
)

link_directories( 
${GAT_SOURCE_DIR}/lib
${GAT_SOURCE_DIR}/external/install/lib
${Boost_LIBRARY_DIRS}
${TBB_LIBRARY_DIRS}
${LAPACK_LIBRARY_DIR}
${BLAS_LIBRARY_DIR}
)

if(HAVE_LOGGER)
    set (LOGGING_LIBS lib_g2logger lib_activeobject)
else()
    set (LOGGING_LIBS )
endif()
#dd_library(pca SHARED ${PCA_LIB_SRCS})

message("TBB_LIBRARIES = ${TBB_LIBRARIES}")
message("Boost_LIBRARIES = ${Boost_LIBRARIES}")

# Build the Sailfish library
add_library(sailfish_core SHARED ${SAILFISH_LIB_SRCS} )
# Build the Sailfish executable
add_executable(sailfish ${SAILFISH_MAIN_SRCS})

# Build the salmon executable
add_executable(salmon ${SALMON_MAIN_SRCS})

add_executable(salmon-read ${SALMON_READ_SRCS})

set_target_properties(salmon-read PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_LIBPTHREAD -D_PBGZF_USE -fopenmp"
    LINK_FLAGS "-DHAVE_LIBPTHREAD -D_PBGZF_USE -fopenmp")

# Link the executable
target_link_libraries(sailfish 
    sailfish_core 
    gff
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    cmph # perfect hashing library
    jellyfish-2.0 
    pthread 
    gomp 
    m
    ${TBB_LIBRARIES}
    ${LOGGING_LIBS}
    ${LIBSAILFISH_LINKER_FLAGS} shark
)

# Link the executable
target_link_libraries(salmon
    sailfish_core 
    gff
    bwa
    clasp
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    cmph # perfect hashing library
    jellyfish-2.0 
    pthread 
    gomp 
    m
    ${TBB_LIBRARIES}
    ${LOGGING_LIBS}
    ${LIBSAILFISH_LINKER_FLAGS} shark
)

add_dependencies(salmon libbwa)

# Link the executable
target_link_libraries(salmon-read
    gff
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARY} 
    cmph # perfect hashing library
    jellyfish-2.0 
    ${CMAKE_SOURCE_DIR}/external/install/lib/libbam.a
    z
    pthread 
    gomp 
    m
    ${TBB_LIBRARIES}
    ${LOGGING_LIBS}
    ${LIBSAILFISH_LINKER_FLAGS} shark
)

##### ======================================

#add_executable(build_transcript_map ${BUILD_TRANSCRIPT_MAP_SRCS})

#target_link_libraries(build_transcript_map
  #sailfish_core
  #${Boost_LIBRARIES}
  #${ZLIB_LIBRARY}
  #cmph
  #jellyfish-2.0
  #pthread
  #gomp
  #m
  #${TBB_LIBRARIES}
  #${LOGGING_LIBS}
  #${LIBSAILFISH_LINKER_FLAGS} shark
#)
   


#add_executable(read_hit_tester ReadHitTester.cpp)
#target_link_libraries(read_hit_tester
    #sailfish_core 
    #gff
    #${Boost_LIBRARIES} 
    #${ZLIB_LIBRARY} 
    #cmph # perfect hashing library
    #jellyfish-2.0 
    #pthread 
    #gomp 
    #m
    #${TBB_LIBRARIES}
    #${LOGGING_LIBS}
    #${LIBSAILFISH_LINKER_FLAGS} shark
#)

#add_executable(BuildKmerLocationIndex BuildKmerLocationIndex.cpp)
#target_link_libraries(BuildKmerLocationIndex
    #sailfish_core 
    #gff
    #${Boost_LIBRARIES} 
    #${ZLIB_LIBRARY} 
    #cmph # perfect hashing library
    #jellyfish-2.0 
    #pthread 
    #gomp 
    #m
    #${TBB_LIBRARIES}
    #${LOGGING_LIBS}
    #${LIBSAILFISH_LINKER_FLAGS} shark
#)


#
#add_executable(performBiasCorrection BiasCorrectionDriver.cpp)
#target_link_libraries(performBiasCorrection
#sailfish_core
#${Boost_LIBRARIES}
#${ZLIB_LIBRARY}
#cmph
#jellyfish-2.0
#pthread
#gomp
#m
#${TBB_LIBRARIES}
#shark
#)


#add_executable(analyze_kmer_tradeoff AnalyzeKmerTradeoff.cpp)
#target_link_libraries(analyze_kmer_tradeoff 
    #sailfish_core 
    #${Boost_LIBRARIES} 
    #${ZLIB_LIBRARY} 
    #cmph # perfect hashing library
    #jellyfish-2.0 
    #pthread 
    #gomp 
    #m
    #${LOGGING_LIBS}
    #${TBB_LIBRARIES}
    #shark
#)

# add_executable(compute_transcript_features ComputeBiasFeatures.cpp)
# target_link_libraries(compute_transcript_features
# ${Boost_LIBRARIES}
# ${ZLIB_LIBRARY}
# jellyfish-2.0
# pthread
# gomp
# m
# ${TBB_LIBRARIES}
# )

# add_executable(perform_bias_correction PerformBiasCorrection.cpp cokus.cpp PCA.cpp)
# target_link_libraries(perform_bias_correction
# ${Boost_LIBRARIES}
# ${ZLIB_LIBRARY}
# pthread
# gomp
# m
# shark # The machine learning library
# )

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${GAT_SOURCE_DIR}" CACHE PATH "Default install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR lib )
set(INSTALL_BIN_DIR bin )
set(INSTALL_INCLUDE_DIR include )

install(DIRECTORY 
        ${GAT_SOURCE_DIR}/external/install/lib/ 
        DESTINATION ${INSTALL_LIB_DIR}
    )

# install(FILES ${Boost_LIBRARIES}
# 	           DESTINATION ${INSTALL_LIB_DIR})

install(TARGETS sailfish sailfish_core
                RUNTIME DESTINATION bin 
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
        )

set(POST_INSTALL_SCRIPT ${GAT_SOURCE_DIR}/cmake/PostInstall.cmake)

install(
    CODE 
    "
    execute_process(COMMAND \"${CMAKE_COMMAND}\"
                            -DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}
                            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} 
                            -P \"${POST_INSTALL_SCRIPT}\")
    "
)

# install(FILES ${GAT_SOURCE_DIR}/scripts/SFPipeline.py DESTINATION scripts/SFPipeline.py )
# install(FILES ${GAT_SOURCE_DIR}/experimental_configs/SRX016368_25mers.cfg DESTINATION experimental_configs/SRX016368_25mers.cfg )

include(InstallRequiredSystemLibraries)
add_test( NAME simple_test COMMAND ${CMAKE_COMMAND} -DTOPLEVEL_DIR=${GAT_SOURCE_DIR} -P ${GAT_SOURCE_DIR}/cmake/SimpleTest.cmake )

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
# 	jellyfish-2.0 pthread 
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

