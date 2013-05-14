set( CMAKE_C_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-Wno-write-strings" )
include_directories( . )

add_library( mlrmcl STATIC coarsen.c match.c ccgraph.c memory.c util.c 
metrics.c merge.c graph.c io.c subgraphs.c myqsort.c mlmcl.c mclbase.c 
mclutils.c rmcl.c)
