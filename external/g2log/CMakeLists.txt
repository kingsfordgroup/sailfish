# CMakeLists.txt cmake configuration for g2log test
# g2log is a KjellKod Logger
# 2011 @author Kjell Hedström, hedstrom@kjellkod.cc
# ==================================================================
# 2010 by KjellKod.cc. This is PUBLIC DOMAIN to use at your own
#    risk and comes  with no warranties.
#
# This code is yours to share, use and modify with no strings attached
#   and no restrictions or obligations.
# ===================================================================

# Below are details for compiling on Windows and Linux
# by default only an example g2log binary is created
# the performance and unit tests creation can be enabled by switching their
# OPTIONs from OFF to ON --- See below at around line 110

# 2012-05-29: justthread is no longer necessary on Windows (vs2011) and
#             linux (gcc4.7.1)


# WINDOWS == README: Example how to setup environment + running an example
# Below written for VS11 (2012)
# 1. please use the "Visual Studio Command Prompt 11 (2012)"
# 2. from the g2log folder
#    mkdir build
#    cd build;
# 3. cmake -DCMAKE_BUILD_TYPE=Release -G "Visual Studio 11" ..
#              the "Visual Studio 11" .. does not require just::thread!
# 4. msbuild g2log_by_kjellkod.sln /p:Configuration=Release
# 5. Release\g2log-FATAL-example.exe
#
#
# . LINUX:To try this out from folder g2log:
#    mkdir build
#    cd build
#    cmake ..    # create makefiles in g2log/build directory
#    make    # link active_object, g2log and example code to get an "example" executable
#     ./g2log-FATAL-example
# ============================================================================

set(LOG_SRC src)

MESSAGE("  LOG_SRC = : ${LOG_SRC}")
include_directories(${LOG_SRC})
SET(ACTIVE_CPP0xx_DIR "Release")

# Detect 64 or 32 bit
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
   # 64-bit project
   SET(64_BIT_OS TRUE)
     MESSAGE("A 64-bit OS detected")
 else()
   SET(64_BIT_OS FALSE)
   MESSAGE("A 32-bit OS detected")
endif()


IF(UNIX)
   MESSAGE("")
   MESSAGE("cmake for *NIX ")
   MESSAGE("if cmake finishes OK, do make")
   MESSAGE("then run './g2log-FATAL-example' or whatever performance test you feel like trying")
   MESSAGE("")
   set(PLATFORM_LINK_LIBRIES rt)
   # justthread not needed from gcc4.7.1
   #           set(PLATFORM_LINK_LIBRIES justthread rt)

   # -rdynamic is needed for correct stack dumps with demangling
   # -D_GLIBCXX_USE_NANOSLEEP is needed for this_thread sleep (unit testing)
   set(CMAKE_CXX_FLAGS "-Wall -rdynamic -Wunused -std=c++11 ${CMAKE_CXX_FLAGS_DEBUG} -pthread -D_GLIBCXX_USE_NANOSLEEP")
   set(SRC_PLATFORM_SPECIFIC ${LOG_SRC}/crashhandler_unix.cpp)

   # -- justthread not needed for gcc4.7.1
   #     set(CMAKE_CXX_FLAGS "-Wall -rdynamic -Wunused -std=c++0x ${CMAKE_CXX_FLAGS_DEBUG} -pthread -I/usr/include/justthread")
   # -- justhtread not needed for gcc4.7.1
   #               include_directories("/usr/include/justthread")
   ENDIF(UNIX)

   # GENERIC STEPS
   # CODE SOURCES these +
   set(SRC_CPP ${LOG_SRC}/g2logworker.cpp ${LOG_SRC}/g2log.cpp ${LOG_SRC}/g2time.cpp)
   set(SRC_H ${LOG_SRC}/g2logworker.h ${LOG_SRC}/g2log.h ${LOG_SRC}/crashhandler.h ${LOG_SRC}/g2time.h ${LOG_SRC}/g2future.h)
   set(SRC_FILES ${SRC_CPP} ${SRC_H} ${SRC_PLATFORM_SPECIFIC})

   # add a ActiveObject library
   add_library(lib_activeobject ${LOG_SRC}/active.cpp ${LOG_SRC}/active.h ${LOG_SRC}/shared_queue.h)
   set_target_properties(lib_activeobject PROPERTIES LINKER_LANGUAGE CXX)

   # add a g2log library
   include_directories(src)
   include_directories(${LOG_SRC})
   MESSAGE("  LOG_SRC = : ${LOG_SRC}")
   MESSAGE("  g2logger files: [${SRC_FILES}]")
   add_library(lib_g2logger  ${SRC_FILES})
   set_target_properties(lib_g2logger PROPERTIES LINKER_LANGUAGE CXX)
   target_link_libraries(lib_g2logger lib_activeobject)


