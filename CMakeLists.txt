cmake_minimum_required (VERSION 2.8)

enable_testing()

project (Sailfish)

set(CPACK_PACKAGE_VERSION "0.6.3")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "6")
set(CPACK_PACKAGE_VERSION_PATCH "3")
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "Carnegie Mellon University")
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
set (CMAKE_CXX_FLAGS "-g -funroll-loops -fPIC -fomit-frame-pointer -Ofast -DHAVE_ANSI_TERM -DHAVE_SSTREAM -DHAVE_CONFIG_H -Wall -std=c++11 -Wreturn-type -Werror=return-type")

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
    #set(HAVE_LOGGER FALSE)
else()
    # We're in sane linux world
   set (SHARED_LIB_EXTENSION "so")
    set (LIBSAILFISH_LINKER_FLAGS "")
endif()

set( BOOST_EXTRA_FLAGS "--layout=tagged" )
##
# Compiler-specific C++11 activation.
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
##
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()

    if (NOT APPLE OR 
            GCC_VERSION VERSION_GREATER 4.8.2 OR
            GCC_VERSION VERSION_EQUAL 4.8.2) 
        #Perhaps the logger now works in both worlds
        set (CMAKE_CXX_FLAGS "-DHAVE_LOGGER ${CMAKE_CXX_FLAGS}")
        set (HAVE_LOGGER TRUE)
    else()
        message(WARNING "Logging requires gcc >= 4.8.2 on OSX; it will be disabled")
    endif()

    set(WARNING_IGNORE_FLAGS "${WARNING_IGNORE_FLAGS} -Wno-unused-local-typedefs")
    set (BOOST_CONFIGURE_TOOLSET "")
    set (BOOST_EXTRA_FLAGS "")
    # Tentatively, we support clang now
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # If we're on OSX using clang, try to use libc++ instead of libstdc++
    if (APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set (BOOST_CONFIGURE_TOOLSET "--with-toolset=clang")
        set (BOOST_EXTRA_FLAGS toolset=clang cxxflags="-stdlib=libc++"  linkflags="-stdlib=libc++")
    endif()
    # There's currently a bug with clang-3.4 & Boost 1.55 -- this hack fixes it
    # but we should do something better (does this break things if CPU doesn't
    # have 128-bit support)?
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_HAS_INT128") 
    # We should have the logger under clang
    set (CMAKE_CXX_FLAGS "-DHAVE_LOGGER ${CMAKE_CXX_FLAGS}")
    set (HAVE_LOGGER TRUE)
 
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

##
#  Update the CXX flags according to the system and compiler
##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_IGNORE_FLAGS}")

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
	message(FATAL_ERROR "zlib must be installed before configuration & building can proceed")
endif()

##
# Set the latest version and look for what we need
##
set(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0" "1.54" "1.55" "1.56")
find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer chrono program_options serialization)
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
    find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer chrono program_options serialization REQUIRED)
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
    set(BOOST_WILL_RECONFIGURE TRUE)
    set(FETCH_BOOST FALSE)
    message("Build system will fetch and build Boost")
    message("==================================================================")
    ExternalProject_Add(libboost
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        URL http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_56_0
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        PATCH_COMMAND patch -p2 < ${CMAKE_CURRENT_SOURCE_DIR}/external/boost156.patch
        CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_56_0/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} ${BOOST_BUILD_LIBS} --prefix=<INSTALL_DIR>
        BUILD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_56_0/b2 -d0 -j4 ${BOOST_EXTRA_FLAGS} cxxflags=${BOOST_CXX_FLAGS} install
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

message("Build system will fetch and build the Cereal serialization library")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libcereal
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/cereal-v1.0.0.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.0.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.0.0/build
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include && cp -r <SOURCE_DIR>/include/cereal <INSTALL_DIR>/include"
)

message("Build system will fetch and build BWA (for Salmon)")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libbwa
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/bwa-master.tar.gz 
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/bwa-master
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND sh -c "make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include/bwa && cp libbwa.a <INSTALL_DIR>/lib && cp *.h <INSTALL_DIR>/include/bwa && cp is.c bwtindex.c bwt_gen.c QSufSort.c ${CMAKE_CURRENT_SOURCE_DIR}/src/"
    BUILD_IN_SOURCE TRUE
)


message("Build system will fetch and build CLASP (for Salmon)")
message("==================================================================")
ExternalProject_Add(libclasp
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/clasp_v1_1.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/clasp_v1_1
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND sh -c "make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} && ar cr libclasp.a <SOURCE_DIR>/libs/*.o"
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/lib && mkdir -p <INSTALL_DIR>/include/clasp && cp libclasp.a <INSTALL_DIR>/lib && cp libs/*.h <INSTALL_DIR>/include/clasp"
    BUILD_IN_SOURCE TRUE
)



message("Build system will fetch and build CMPH")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libcmph
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/cmph-2.0.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cmph-2.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/cmph-2.0/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND ${MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} 
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
)

message("Build system will fetch and build Jellyfish")
message("==================================================================")
ExternalProject_Add(libjellyfish
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL ftp://ftp.genome.umd.edu/pub/jellyfish/jellyfish-2.1.3.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.1.3
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.1.3/configure --prefix=<INSTALL_DIR> CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND ${MAKE} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} 
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install && 
                    cp config.h <INSTALL_DIR>/include/jellyfish-2.1.3/jellyfish/ &&
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
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(TBB_COMPILER "clang")
else()
    set(TBB_COMPILER "gcc")
endif()

message("Build system will fetch and build Intel Threading Building Blocks")
message("==================================================================")
# These are useful for the custom install step we'll do later
set(TBB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb42_20140601)
set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
ExternalProject_Add(libtbb
	DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb42_20140601oss_src.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb42_20140601
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND "${TBB_PATCH_STEP}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS 
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
# If we're fetching tbb, we need to have dummy paths for these variables
# so that CMake won't complain
##
if(TBB_WILL_RECONFIGURE)
    set(TBB_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(TBB_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
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

###
# Shark library
####
if (NOT Boost_FOUND)
    # Shark depends on Boost.  If it wasn't already found, then 
    # add a dependency between the Boost we fetch and this target.
    add_dependencies(libshark libboost)
endif()

message("Build system will fetch and build Shark machine learning library")
message("==================================================================")
# These are useful for the custom install step we'll do later
set(SHARK_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Shark)
set(SHARK_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
ExternalProject_Add(libshark
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/Shark.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Shark
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CMAKE_ARGS -DOPT_DYNAMIC_LIBRARY=FALSE -DBoost_NO_SYSTEM_PATHS=TRUE -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS} 
               -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS} -DOPT_MAKE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${SHARK_INSTALL_DIR}
               -DCMAKE_CXX_FLAGS=-std=c++11 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    BUILD_COMMAND sh -c "make"
    INSTALL_COMMAND sh -c "make install"
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
)

if(HAVE_LOGGER)
message("Build system will compile g2log logging library")
message("==================================================================")
ExternalProject_Add(libg2log
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
  URL http://www.cs.cmu.edu/~robp/files/g2log.tar.gz
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log
  INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
  BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log/build
  #UPDATE_COMMAND sh -c "mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log/g2log/build"
  UPDATE_COMMAND sh -c "mv ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log/g2log/* ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log && mkdir ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log/build"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/external/install #-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  BUILD_COMMAND sh -c "make"
  INSTALL_COMMAND sh -c "cp <SOURCE_DIR>/src/*.h <INSTALL_DIR>/include/ && cp <BINARY_DIR>/*.a <BINARY_DIR>/*.${SHARED_LIB_EXTENSION} <INSTALL_DIR>/lib/"
)
endif()


message("Build system will compile libgff")
message("==================================================================")
ExternalProject_Add(libgff
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/libgff.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libgff/build
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}/external/install
)

#message("Build system will compile parallel samtools")
#message("==================================================================")
#ExternalProject_Add(libsamtools
    #DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    #URL http://www.cs.cmu.edu/~robp/files/samtools.tgz
    #SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
    #BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
    #INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    #CONFIGURE_COMMAND ""
    #BUILD_COMMAND $(MAKE)
    #INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include/samtools && mkdir -p <INSTALL_DIR>/include/samtools/pbgzip && mkdir -p <INSTALL_DIR>/include/samtools/bcftools && cp <SOURCE_DIR>/*.h <INSTALL_DIR>/include/samtools/ &&  cp <SOURCE_DIR>/pbgzip/*.h <INSTALL_DIR>/include/samtools/pbgzip/ && cp <SOURCE_DIR>/bcftools/*.h <INSTALL_DIR>/include/samtools/bcftools/ &&  cp <SOURCE_DIR>/libbam.a <INSTALL_DIR>/lib/"
#)

message("Build system will compile htslib")
message("==================================================================")
ExternalProject_Add(libhts
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/samtools/htslib/archive/1.1.tar.gz -o htslib-1.1.tar.gz &&
                     tar -xzf htslib-1.1.tar.gz &&
                     rm -fr htslib &&
                     mv -f htslib-1.1 htslib
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/htslib
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install prefix=<INSTALL_DIR> 
)

message("Build system will compile samtools")
message("==================================================================")
ExternalProject_Add(libsamtools
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    DOWNLOAD_COMMAND curl -k -L https://github.com/samtools/samtools/archive/1.1.tar.gz -o samtools-1.1.tar.gz &&
                     tar -xzf samtools-1.1.tar.gz &&
                     rm -fr samtools &&
                     mv -f samtools-1.1 samtools
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
    CONFIGURE_COMMAND ""
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    BUILD_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install prefix=<INSTALL_DIR> &&
                    mv <SOURCE_DIR>/libbam.a <INSTALL_DIR>/lib/libbam.a
)


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

