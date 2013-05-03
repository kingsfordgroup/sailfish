cmake_minimum_required (VERSION 2.6)

project (Sailfish)

set (CMAKE_CXX_COMPILER g++)

set (CMAKE_CXX_FLAGS "-O3 -DHAVE_ANSI_TERM -DHAVE_SSTREAM -DHAVE_CONFIG_H -Wall -std=c++11")

set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set (CMAKE_VERBOSE_MAKEFILE true)

if ( DEFINED CUSTOM_BOOST_PATH )
	set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CUSTOM_BOOST_PATH})
	set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CUSTOM_BOOST_PATH}/lib )
endif ( DEFINED CUSTOM_BOOST_PATH )

set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)

set(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0")
find_package(Boost 1.53.0 COMPONENTS iostreams filesystem system thread timer program_options)

if (NOT Boost_FOUND)
	message(FATAL_ERROR 
		"Sailfish cannot be compiled without Boost.\nPlease visit http://www.boost.org/ and install Boost")
endif()


find_package (ZLIB)

###
#
# Build external dependencies (except Boost)
#
###

set(EXTERNAL_LIBRARY_PATH $CMAKE_CURRENT_SOURCE_DIR/lib)

message(status "Building and Installing CMPH")
message(status "============================")
include(ExternalProject)
ExternalProject_Add(libcmph
	DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/cmph-2.0
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/cmph-2.0/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
)
message(status "CMPH installed.")

message(status "Downloading, Building and Installing Jellyfish")
message(status "============================")
ExternalProject_Add(libjellyfish
	DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
	URL http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.10.tar.gz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-1.1.10
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/external/jellyfish-1.1.10/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install && cp config.h <INSTALL_DIR>/include/jellyfish-1.1.10/
)
message(status "Jellyfish installed.")

message(status "Downloading, Building and Installing Intel TBB")
message(status "============================")

# These are useful for the custom install step we'll do later
set(TBB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss)
set(TBB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install)

ExternalProject_Add(libtbb
	DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external
	URL http://threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb41_20130314oss_src.tgz
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/tbb41_20130314oss
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make cfg=release tbb_build_prefix=LIBS
    INSTALL_COMMAND sh -c "cp ${TBB_SOURCE_DIR}/build/LIBS_release/*.so* ${TBB_INSTALL_DIR}/lib && cp -r ${TBB_SOURCE_DIR}/include/* ${TBB_INSTALL_DIR}/include"
    BUILD_IN_SOURCE 1
)
message(status "TBB installed.")

add_subdirectory ( external )

###
#
# Done building external dependencies.
#
###

# Recurse into Sailfish source directory
add_subdirectory ( src )
