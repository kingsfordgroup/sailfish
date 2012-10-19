set( CMAKE_C_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-fPIC")
include_directories( . )

add_library( nnls SHARED nnls.c )
