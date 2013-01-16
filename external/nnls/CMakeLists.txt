set( CMAKE_C_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-fPIC -")
include_directories( include )

#add_library( nnls SHARED nnls.cc )
add_library( nnls SHARED G12.c  H12.c  myMath.c  nnls.c) 
