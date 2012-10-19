cmake_minimum_required (VERSION 2.6)

project (Sailfish)

set (CMAKE_CXX_COMPILER g++)#-mp-4.7)

set (CMAKE_CXX_FLAGS "-g -O3 -DHAVE_JFCONFIG_H -std=c++11")

set (GAT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set (CMAKE_VERBOSE_MAKEFILE true)

if ( DEFINED CUSTOM_BOOST_PATH )
	set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CUSTOM_BOOST_PATH})
	set (CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CUSTOM_BOOST_PATH}/lib )
endif ( DEFINED CUSTOM_BOOST_PATH )

set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.50.0 COMPONENTS iostreams filesystem system )
find_package (ZLIB)

if ( DEFINED SEQAN_PATH )
   message ("Using SeqAn from ${SEQAN_PATH}")
else ( DEFINED SEQAN_PATH )
     message ("You must provide a path for SeqAn (e.g. cmake .. -DSEQAN_PATH=...)")
endif ( DEFINED SEQAN_PATH )

add_subdirectory ( external )
add_subdirectory ( src )
