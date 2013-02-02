set( CMAKE_C_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-Wno-write-strings -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -fopenmp")
include_directories( . )

add_library( shotgun STATIC 
lasso.cpp
logreg.cpp
common.h 
mmio.c 
mmio.h 
read_matrix_market.cpp 
write_matrix_market.cpp 
shared.cpp)
