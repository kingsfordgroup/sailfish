cmake_minimum_required (VERSION 2.8)
enable_testing()

project (Sailfish)

set(CPACK_PACKAGE_VERSION "0.10.0")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "10")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "Stony Brook University / Carnegie Mellon University")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Sailfish - Alignment-free RNA-seq isoform quantification")
set(CPACK_PACKAGE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-Source")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set (WARNING_IGNORE_FLAGS "-Wno-deprecated-register")
set (BOOST_CXX_FLAGS "-Wno-deprecated-register -std=c++11")
## Prefer static to dynamic libraries
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

## Set the standard required compile flags
# Nov 18th --- removed -DHAVE_CONFIG_H
set (CMAKE_CXX_FLAGS "-pthread -funroll-loops -fPIC -fomit-frame-pointer -Ofast -DHAVE_ANSI_TERM -DHAVE_SSTREAM -DRAPMAP_SALMON_SUPPORT -Wall -std=c++11 -Wreturn-type -Werror=return-type")

##
# OSX is strange (some might say, stupid in this regard).  Deal with it's quirkines here.
##
if (APPLE)
    # To allow ourselves to build a dynamic library, we have to tell the compiler
    # that, yes, the symbols will be around at runtime.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
    set (LIBSAILFISH_LINKER_FLAGS "-all_load")
    # In order to "think different", we also have to use non-standard suffixes
    # for our shared libraries
    set(SHARED_LIB_EXTENSION "dylib")
else()
    # We're in sane linux world
   set (SHARED_LIB_EXTENSION "so")
   set (LIBSAILFISH_LINKER_FLAGS "")
endif()

set( BOOST_EXTRA_FLAGS "--layout=tagged" )
## this get's set differently below if we
## are on clang & apple
set (NON_APPLECLANG_LIBS gomp rt)
set (PTHREAD_LIB)

##
# Compiler-specific C++11 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
##
# First take care of what to do if we have gcc
##
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    # If we're on OSX
    if (APPLE AND NOT (GCC_VERSION VERSION_GREATER 4.8.2 OR GCC_VERSION VERSION_EQUAL 4.8.2))
        message(FATAL_ERROR "When building under OSX, ${PROJECT_NAME} requires "
                            "either clang or g++ >= 4.8.2")
    elseif (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()

    set (GCC TRUE)

    # Put complete static linking on hold for the time-being
    # If we're not on OSX, make an attempt to compile everything statically
    #if (NOT APPLE)
    #set (CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
    #set (CMAKE_EXE_LINK_FLAGS "-static")
    set (PTHREAD_LIB "pthread")
    #endif()

    # If we're on Linux (i.e. not OSX) and we're using
    # gcc, then set the -static-libstdc++ flag
    if (NOT APPLE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
    endif()

    set (WARNING_IGNORE_FLAGS "${WARNING_IGNORE_FLAGS} -Wno-unused-local-typedefs")
    set (BOOST_TOOLSET "gcc")
    set (BOOST_CONFIGURE_TOOLSET "--with-toolset=gcc")
	set (BCXX_FLAGS "-std=c++11")
    set (BOOST_EXTRA_FLAGS toolset=gcc cxxflags=${BCXX_FLAGS})
# Tentatively, we support clang now
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CLANG TRUE)
    # If we have libc++, then try and use it
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-stdlib=libc++ HAVE_LIBCPP)
    if (HAVE_LIBCPP)
        message ("It appears that you're compiling with clang and that libc++ is available, so I'll use that")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	    set (BOOST_TOOLSET "clang")
        set (BOOST_CONFIGURE_TOOLSET "--with-toolset=clang")
	    set (BCXX_FLAGS "-stdlib=libc++ -DBOOST_HAS_INT128")
	    set (BOOST_EXTRA_FLAGS toolset=clang cxxflags=${BCXX_FLAGS} linkflags="-stdlib=libc++")
        set (JELLYFISH_CXX_FLAGS "-stdlib=libc++")
    # Otherwise, use libstdc++ (and make it static)
    else()
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
    endif()
    # There's currently a bug with clang-3.4 & Boost 1.55 -- this hack fixes it
    # but we should do something better (does this break things if CPU doesn't
    # have 128-bit support)?
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_HAS_INT128")

    if (APPLE)
        set (NON_APPLECLANG_LIBS "")
    else()
        set (PTHREAD_LIB "pthread")
    endif()
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

## TODO: Figure out how to detect this automatically
# If the "assembler" is too old, tell TBB not to compile
# with -mrtm
if (NO_RTM)
    set (TBB_CXXFLAGS "-mno-rtm")
endif()

include(ExternalProject)

##
#  Update the CXX flags according to the system and compiler
##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_IGNORE_FLAGS}")

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message ("Making Debug build")
    set (CMAKE_CXX_FLAGS_DEBUG "-g ${CMAKE_CXX_FLAGS}")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
    message ("Making Release build")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS}")
else ()
    message ("Making Default build type")
endif ()

##
# Record this top-level path
##
set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Have CMake tell us what it's doing
set (CMAKE_VERBOSE_MAKEFILE true)

##
# Super-secret override
##
if ( DEFINED CUSTOM_BOOST_PATH )
	set (CMAKE_INCLUDE_PATH ${CUSTOM_BOOST_PATH} ${CMAKE_INCLUDE_PATH})
    set (CMAKE_LIBRARY_PATH ${CUSTOM_BOOST_PATH}/lib ${CMAKE_LIBRARY_PATH})
endif ( DEFINED CUSTOM_BOOST_PATH )

##
# We want static, multithreaded boost libraries
##
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED ON)
#set (Boost_USE_STATIC_RUNTIME OFF)

find_package (ZLIB)
if (NOT ZLIB_FOUND)
	message (FATAL_ERROR "zlib must be installed before configuration & building can proceed")
endif()

##
# Set the latest version and look for what we need
##
set(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0" "1.54" "1.55" "1.56" "1.57.0" "1.58" "1.59")
find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer chrono program_options serialization)
message("BOOST_INCLUDEDIR = ${BOOST_INCLUDEDIR}")
message("BOOST_LIBRARYDIR = ${BOOST_LIBRARYDIR}")
message("Boost_FOUND = ${Boost_FOUND}")
include(ExternalProject)

##
#  If we had to fetch Boost, the reconfigure step will re-run cmake.  The second configuration
#  pass is executed with the BOOST_RECONFIGURE flag set. This should allow our newly
#  installed Boost to be found by CMake.
##
if (BOOST_RECONFIGURE)
    message("Executing Boost Reconfiguration")
    unset(Boost_FOUND CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    unset(Boost_LIBRARY_DIRS CACHE)
    unset(Boost_LIBRARIES CACHE)
    unset(BOOST_ROOT CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)

    set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer chrono program_options REQUIRED)
    set(FETCH_BOOST FALSE)
endif()



##
# Either inform the user of how to obtain Boost, or, if they passed in the FETCH_BOOST
# option, go and grab it for them.
##
if ((NOT Boost_FOUND) AND (NOT FETCH_BOOST))
	message(FATAL_ERROR
		"Sailfish cannot be compiled without Boost.\n"
        "It is recommended to visit http://www.boost.org/ and install Boost according to those instructions.\n"
        "This build system can also download and install a local version of boost for you (this takes a lot of time).\n"
        "To fetch and build boost locally, call cmake with -DFETCH_BOOST=TRUE"
    )
elseif(FETCH_BOOST)
    ## Let the rest of the build process know we're going to be fetching boost
    set (BOOST_LIB_SUBSET --with-atomic --with-chrono --with-container --with-date_time --with-exception
                          --with-filesystem --with-math
                          --with-program_options --with-system --with-thread --with-iostreams
                          --with-timer)
    set (BOOST_WILL_RECONFIGURE TRUE)
    set (FETCH_BOOST FALSE)
    message("Build system will fetch and build Boost")
    message("==================================================================")
    ExternalProject_Add(libboost
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        DOWNLOAD_COMMAND curl -k -L  http://downloads.sourceforge.net/project/boost/boost/1.59.0/boost_1_59_0.tar.gz -o boost_1_59_0.tar.gz &&
                               tar xzf boost_1_59_0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_59_0
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        #PATCH_COMMAND patch -p2 < ${CMAKE_CURRENT_SOURCE_DIR}/external/boost156.patch
	CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_59_0/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} ${BOOST_BUILD_LIBS} --prefix=<INSTALL_DIR>
	BUILD_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_59_0/b2 -d0 -j2 ${BOOST_LIB_SUBSET} toolset=${BOOST_TOOLSET} ${BOOST_EXTRA_FLAGS} cxxflags=${BOOST_CXX_FLAGS} link=static install
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
    )

    ##
    # After we've installed boost,
    ##
    SET( RECONFIG_FLAGS ${RECONFIG_FLAGS} -DBOOST_WILL_RECONFIGURE=FALSE -DBOOST_RECONFIGURE=TRUE -DFETCH_BOOST=FALSE)
    ExternalProject_Add_Step(libboost reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
        DEPENDEES install
    )
endif()

##
# If we're fetching boost and we need to have dummy paths for these variables
# so that CMake won't complain
##
if (BOOST_WILL_RECONFIGURE)
    message("Setting Temporary Boost paths")
    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    set(Boost_FOUND TRUE)
endif()


message("BOOST INCLUDE DIR = ${Boost_INCLUDE_DIR}")
message("BOOST INCLUDE DIRS = ${Boost_INCLUDE_DIRS}")
message("BOOST LIB DIR = ${Boost_LIBRARY_DIRS}")
message("BOOST LIBRAREIS = ${Boost_LIBRARIES}")

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

message("Build system will build libdivsufsort")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libdivsufsort
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort.zip
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort-master
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libdivsufsort-master/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DBUILD_DIVSUFSORT64=TRUE -DUSE_OPENMP=TRUE -DBUILD_SHARED_LIBS=FALSE
    )
set(SUFFARRAY_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)


message("Build system will fetch and build Jellyfish")
message("==================================================================")
ExternalProject_Add(libjellyfish
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz -o jellyfish-2.2.3.tgz &&
	rm -fr jellyfish-2.2.3 &&
    	tar -xzvf jellyfish-2.2.3.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.2.3
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.2.3/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${JELLYFISH_CXX_FLAGS}
    BUILD_COMMAND ${MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${JELLYFISH_CXX_FLAGS}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install &&
                    cp config.h <INSTALL_DIR>/include/jellyfish-2.2.3/jellyfish/ &&
                    cp config.h <INSTALL_DIR>/include/
)

find_package(TBB)
##
#
# Fetch and build Intel's Threading Building Blocks library.
#
##
if(NOT TBB_FOUND)

set(TBB_WILL_RECONFIGURE TRUE)
# Set the appropriate compiler
if (CLANG)
    set(TBB_COMPILER "clang")
else()
    set(TBB_COMPILER "gcc")
endif()

message("Build system will fetch and build Intel Threading Building Blocks")
message("==================================================================")
# These are useful for the custom install step we'll do later
set(TBB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb43_20140724oss)
set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)

if ("${TBB_COMPILER}" STREQUAL "gcc")
    ## Don't know why it's a problem yet, but if we're using
    ## GCC, get rid of the DO_ITT_NOTIFY flag
    set(TBB_CXXFLAGS "${TBB_CXXFLAGS} -UDO_ITT_NOTIFY")
endif()

ExternalProject_Add(libtbb
	DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL https://distfiles.macports.org/tbb/tbb43_20140724oss_src.tgz
    DOWNLOAD_COMMAND curl -k -L  https://distfiles.macports.org/tbb/tbb43_20140724oss_src.tgz -o tbb_20140724oss_src.tgz &&
                     tar -xzvf tbb_20140724oss_src.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb43_20140724oss
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND "${TBB_PATCH_STEP}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CXXFLAGS=${TBB_CXXFLAGS} lambdas=1 compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS
    INSTALL_COMMAND sh -c "cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.${SHARED_LIB_EXTENSION}* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
    BUILD_IN_SOURCE 1
)

SET( RECONFIG_FLAGS ${RECONFIG_FLAGS} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE)
ExternalProject_Add_Step(libtbb reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} ${RECONFIG_FLAGS}
        DEPENDEES install
)
endif()

## 
#  Find HDF5, and prefer the static libs
##
#set(HDF5_USE_STATIC_LIBRARIES TRUE)
#find_package(HDF5 REQUIRED)

#if (HDF5_FOUND)
#    include_directories(${HDF5_INCLUDE_DIRS})
#    #set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
#    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed -ldl")
#    #set_target_properties(sailfish sailfish_core PROPERTIES LINK_SEARCH_END_STATIC 1)
#else()
#    message(FATAL_ERROR "HDF5 libraries are required but were not found!")
#endif()
##
# Done with HDF5
##

##
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(TBB_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    #set(TBB_LIBRARIES tbb tbbmalloc)
    set(TBB_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libtbb.${SHARED_LIB_EXTENSION}
                      ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libtbbmalloc.${SHARED_LIB_EXTENSION}
    )
endif()

##
#  Similar to the Boost trick above, the libtbb reconfigure should force this code
#  to be run on the second configuration pass, where it should appropriately set the
#  TBB_INSTALL_DIR variable.
##
if (TBB_RECONFIGURE)
    unset(TBB_FOUND CACHE)
    unset(TBB_INSTALL_DIR CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    message("TBB_INSTALL_DIR = ${TBB_INSTALL_DIR}")
    find_package(TBB)
endif()


message("TBB_LIBRARIES = ${TBB_LIBRARIES}")

message("Build system will compile libgff")
message("==================================================================")
ExternalProject_Add(libgff
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/Kingsford-Group/libgff/archive/v1.0.tar.gz -o libgff.tgz &&
    	tar -xzvf libgff.tgz &&
        rm -fr libgff &&
	    mv libgff-1.0 libgff
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/external/install
)

set (FAST_MALLOC_LIB "")
set (HAVE_FAST_MALLOC FALSE)

# See if we have Jemalloc
find_package(Jemalloc)
if (Jemalloc_FOUND)
    message("Found Jemalloc library --- using this memory allocator")
    set (FAST_MALLOC_LIB ${JEMALLOC_LIBRARIES})
    set (HAVE_FAST_MALLOC TRUE)
endif()

if (NOT HAVE_FAST_MALLOC)
    # See if we have Tcmalloc
    find_package(Tcmalloc)
    if (Tcmalloc_FOUND)
        message("Fount TCMalloc library --- using this memory allocator")
        set (TCMALLOC_LIB ${Tcmalloc_LIBRARIES})
        set (FAST_MALLOC_LIB ${TCMALLOC_LIB})
        set (HAVE_FAST_MALLOC TRUE)
    endif()
endif()

if (NOT HAVE_FAST_MALLOC)
    message("Build system will fetch and use JEMalloc")
    message("==================================================================")
    ExternalProject_Add(libjemalloc
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/jemalloc/archive/3.6.0.tar.gz -o jemalloc-3.6.0.tar.gz &&
        tar -xzf jemalloc-3.6.0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jemalloc-3.6.0
        BUILD_IN_SOURCE TRUE
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        CONFIGURE_COMMAND sh -c "CC=${CMAKE_C_COMPILER} ./autogen.sh --prefix=<INSTALL_DIR>"
        INSTALL_COMMAND cp -r lib <INSTALL_DIR>/ && cp -r include <INSTALL_DIR>/
        )

    set (FAST_MALLOC_LIB ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib/libjemalloc.a)
    set (HAVE_FAST_MALLOC TRUE)
endif ()

##
## This depenency is for RapMap
##
message("Build system will fetch and build SparseHash")
message("==================================================================")
ExternalProject_Add(libsparsehash
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/COMBINE-lab/sparsehash/archive/sparsehash-2.0.2.tar.gz -o sparsehash-2.0.2.tar.gz &&
        tar -xzf sparsehash-2.0.2.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/sparsehash-sparsehash-2.0.2
    BUILD_IN_SOURCE TRUE
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND sh -c "CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ./configure --prefix=<INSTALL_DIR>"
    INSTALL_COMMAND make install
)

###
#
# Done building external dependencies.
#
###

###
#
#  Grab RapMap sources for quasi-mapping code --- DURING CONFIGURE TIME!
#
####
if(NOT FETCHED_RAPMAP)
  exec_program(${CMAKE_CURRENT_SOURCE_DIR}/scripts/fetchRapMap.sh)
  set(FETCHED_RAPMAP TRUE CACHE BOOL "Has RapMap been fetched?" FORCE)
endif()

###
#
# Done building external dependencies.
#
###

set (CPACK_SOURCE_IGNORE_FILES
"/src/PCA.cpp"
"/src/PCAUtils.cpp"
"/build/"
"/scripts/AggregateToGeneLevel.py"
"/scripts/ExpressionTools.py"
"/scripts/GenerateExpressionFiles.sh"
"/scripts/ParseSoftFile.py"
"/scripts/PlotCorrelation.py"
"/scripts/junk"
"/scripts/sfstrace.log"
"/scripts/SFPipeline.py"
"/bin/"
"/lib/"
"/sample_data/"
"PublishREADMEToWebsite.sh"
"/external/"
"/src/obsolete/"
"/include/obsolete/"
"WebsiteHeader.txt"
"/experimental_configs/"
".git/")

message("CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")

# Recurse into Sailfish source directory
add_subdirectory ( src )

# build a CPack driven installer package
include (CPack)
