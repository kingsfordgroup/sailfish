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

## Set the standard required compile flags
set (CMAKE_CXX_FLAGS "-g -funroll-loops -fPIC -O3 -DHAVE_ANSI_TERM -DHAVE_SSTREAM -DHAVE_CONFIG_H -Wall -std=c++11")

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


    set (BOOST_CONFIGURE_TOOLSET "")
    set (BOOST_EXTRA_FLAGS "")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(FATAL_ERROR "Currently, the only supported compiler is GCC(>=4.7).  We hope to support Clang soon.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set (BOOST_CONFIGURE_TOOLSET "--with-toolset=clang")
    set (BOOST_EXTRA_FLAGS toolset=clang cxxflags="-stdlib=libc++"  linkflags="-stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
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
# We want dyanmic, multithreaded boost libraries
##
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)

# unset(Boost_FOUND CACHE)
# unset(Boost_INCLUDE_DIR CACHE)
# unset(Boost_INCLUDE_DIRS CACHE)
# unset(Boost_LIBRARY_DIRS CACHE) 
# message("BOOST_ROOT = ${BOOST_ROOT}")

##
# Set the latest version and look for what we need
##
set(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0")
find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer program_options serialization)

include(ExternalProject)

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
    message("Build system will fetch and build Boost")
    message("==================================================================")
    ExternalProject_Add(libboost
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
        URL http://www.cs.cmu.edu/~robp/files/boost-1.53.0.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_53_0
        INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
        CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_53_0/bootstrap.sh ${BOOST_CONFIGURE_TOOLSET} --prefix=<INSTALL_DIR>
        BUILD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_53_0/b2 ${BOOST_EXTRA_FLAGS} install
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
    )

    ##
    # After we've installed boost, 
    ##
    ExternalProject_Add_Step(libboost reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} -DBOOST_WILL_RECONFUGRE=FALSE -DBOOST_RECONFIGURE=TRUE 
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

##
#  If we had to fetch Boost, the reconfigure step will re-run cmake.  The second configuration
#  pass is executed with the BOOST_RECONFIGURE flag set. This should allow our newly
#  installed Boost to be found by CMake.
##
if (BOOST_RECONFIGURE)
    unset(Boost_FOUND CACHE)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_INCLUDE_DIRS CACHE)
    unset(Boost_LIBRARY_DIRS CACHE) 
    unset(Boost_LIBRARIES CACHE)
    unset(BOOST_ROOT CACHE)
    unset(CMAKE_PREFIX_PATH CACHE)
    set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/include)
    set(Boost_LIBRARY_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/external/install/lib)
    find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer program_options serialization)
endif()


message("BOOST INCLUDE DIR = ${Boost_INCLUDE_DIR}")
message("BOOST INCLUDE DIRS = ${Boost_INCLUDE_DIRS}")
message("BOOST LIB DIR = ${Boost_LIBRARY_DIRS}")
message("BOOST LIBRAREIS = ${Boost_LIBRARIES}")
find_package (ZLIB)

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

##
# Logging library; currently unused so the dependency is temporarily removed
##
# message(status "Building and Installing g2log")
# message(status "============================")
# include(ExternalProject)
# ExternalProject_Add(libg2log
#     DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
#     URL http://www.cs.cmu.edu/~robp/files/g2log.tgz
#     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/g2log
#     INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
#     CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
#     BUILD_COMMAND make
#     #BUILD_IN_SOURCE 1
#     INSTALL_COMMAND sh -c "cp <BINARY_DIR>/*.a <INSTALL_DIR>/lib && cp <SOURCE_DIR>/src/*.h <INSTALL_DIR>/include"
# )
# message(status "g2log installed.")

message("Build system will fetch and build the Cereal serialization library")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libccereal
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL https://github.com/USCiLab/cereal/archive/v1.0.0.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal-1.0.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    UPDATE_COMMAND sh -c "mkdir -p <SOURCE_DIR>/build"
    CMAKE_ARGS -DBOOST_ROOT=/junk -DBoost_NO_SYSTEM_PATHS=TRUE ${INSTALL_DIR} -DBoost_FOUND=FALSE
)

message("Build system will fetch and build CMPH")
message("==================================================================")
include(ExternalProject)
ExternalProject_Add(libcmph
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/cmph-2.0.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cmph-2.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/cmph-2.0/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
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
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-2.1.3/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
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
    ##
    # On OSX, tbb currently has a linker bug so we have to patch up the following file
    ##
   if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    set(TBB_LINK_SYMBOLS "${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss/src/tbb/mac64-tbb-export.lst")
    set(TBB_PATCH_STEP awk "!/__TBB_SYMBOL\( _ZTVN3tbb13tbb_exceptionE \)/" ${TBB_LINK_SYMBOLS} > ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss/tmplnk && mv ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss/tmplnk ${TBB_LINK_SYMBOLS})
   else("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
    set(TBB_PATCH_STEP "")
   endif("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")

else()

    set(TBB_COMPILER "gcc")

endif()

message("Build system will fetch and build Intel Threading Building Blocks")
message("==================================================================")
# These are useful for the custom install step we'll do later
set(TBB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss)
set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
ExternalProject_Add(libtbb
	DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
	URL http://threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb41_20130314oss_src.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    PATCH_COMMAND "${TBB_PATCH_STEP}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make compiler=${TBB_COMPILER} cfg=release tbb_build_prefix=LIBS
    INSTALL_COMMAND sh -c "cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.${SHARED_LIB_EXTENSION}* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
    BUILD_IN_SOURCE 1
)

ExternalProject_Add_Step(libtbb reconfigure
        COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR} -DTBB_WILL_RECONFIGURE=FALSE -DTBB_RECONFIGURE=TRUE
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
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/external/install
  BUILD_COMMAND sh -c "make"
  INSTALL_COMMAND sh -c "cp <SOURCE_DIR>/src/*.h <INSTALL_DIR>/include/ && cp <BINARY_DIR>/*.a <INSTALL_DIR>/lib/"
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

message("Build system will compile parallel samtools")
message("==================================================================")
ExternalProject_Add(libsamtools
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
    URL http://www.cs.cmu.edu/~robp/files/samtools.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/samtools
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND sh -c "mkdir -p <INSTALL_DIR>/include/samtools && mkdir -p <INSTALL_DIR>/include/samtools/pbgzip && mkdir -p <INSTALL_DIR>/include/samtools/bcftools && cp <SOURCE_DIR>/*.h <INSTALL_DIR>/include/samtools/ &&  cp <SOURCE_DIR>/pbgzip/*.h <INSTALL_DIR>/include/samtools/pbgzip/ && cp <SOURCE_DIR>/bcftools/*.h <INSTALL_DIR>/include/samtools/bcftools/ &&  cp <SOURCE_DIR>/libbam.a <INSTALL_DIR>/lib/"
)

# ###
# # Armadillo library
# ####
# if (NOT Boost_FOUND)
#     # Armadillo depends on Boost.  If it wasn't already found, then 
#     # add a dependency between the Boost we fetch and this target.
#     add_dependencies(libarmadillo libboost)
# endif()

# find_package( LAPACK REQUIRED )
# find_package( BLAS REQUIRED )

# message("Build system will fetch and build the Armadillo library")
# message("==================================================================")
# # These are useful for the custom install step we'll do later
# set(ARMADILLO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Armadillo)
# set(ARMADILLO_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)
# ExternalProject_Add(libarmadillo
#     DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
#     URL http://sourceforge.net/projects/arma/files/armadillo-3.900.6.tar.gz
#     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Armadillo
#     INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
#     CMAKE_ARGS -DBoost_NO_SYSTEM_PATHS=TRUE -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIRS} 
#                -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIRS} -DARMA_USE_BLAS=TRUE -DARMA_USE_LAPACK=TRUE
#                -DCMAKE_INSTALL_PREFIX=${ARMADILLO_INSTALL_DIR}
#                -DCMAKE_LIBRARY_PATH=/home/robp/dist/lib

#     BUILD_COMMAND sh -c "make"
#     INSTALL_COMMAND sh -c "make install"
#     BUILD_IN_SOURCE 1
#     UPDATE_COMMAND ""
# )

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

