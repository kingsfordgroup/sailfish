#
# Makefile.am for /tsnnls directory
#

AM_CFLAGS    = -Wall 
SUBDIRS = taucs_basic

EXTRA_DIST = testfull full_example_A.mat full_example_b.mat full_example_x.mat \
	     sp_test1.sh A_01.sparse x_01.mat b_01.mat \
	     sp_test2.sh A_02.sparse x_02.mat b_02.mat \
	     sp_test3.sh A_03.sparse x_03.mat b_03.mat \
	     sp_test4.sh A_04.sparse x_04.mat b_04.mat \
	     sp_test5.sh A_05.sparse x_05.mat b_05.mat \
	     sp_test6.sh A_06.sparse x_06.mat b_06.mat \
	     rrtest.sh rrtestA.sparse rrtestb.mat rrtestx.mat 

pkginclude_HEADERS = tsnnls.h lsqr.h
noinst_HEADERS = tsnnls_blas_wrappers.h taucs_basic/taucs_common.h acint32_type.h
nobase_pkginclude_HEADERS = taucs_basic/taucs.h taucs_basic/colamd.h taucs_basic/f2c.h taucs_basic/taucs_config_build.h taucs_basic/taucs_config_tests.h taucs_basic/taucs_private.h

lib_LTLIBRARIES = libtsnnls.la

libtsnnls_la_SOURCES = taucs_basic/taucs_logging.c \
	taucs_basic/taucs_malloc.c \
	taucs_basic/taucs_ccs_order.c \
	taucs_basic/taucs_ccs_ops.c \
	taucs_basic/taucs_vec_base.c \
	taucs_basic/taucs_complex.c \
	taucs_basic/colamd.c \
	taucs_basic/amdbar.c \
	taucs_basic/amdexa.c \
	taucs_basic/amdtru.c \
	taucs_basic/genmmd.c \
	taucs_basic/taucs_timer.c \
	taucs_basic/taucs_sn_llt.c \
	taucs_basic/taucs_ccs_base.c \
	tlsqr.c \
        tsnnls.c \
	lsqr.c

libtsnnls_la_CFLAGS = $(AM_CFLAGS)

if ARGTABLE 
  bin_PROGRAMS = tsnnls_test
  tsnnls_test_SOURCES = tsnnls_test.c
  tsnnls_test_LDFLAGS = libtsnnls.la
endif

if ARGTABLE
  TESTS_ENVIRONMENT = argtable=1
  SHELL_TESTS = testfull transpose_test sp_test1.sh sp_test2.sh sp_test3.sh sp_test4.sh sp_test5.sh sp_test6.sh rrtest.sh
endif

if HAVE_LAPACK 
  check_PROGRAMS = transpose_test f77blas_test genb_test
  genb_test_SOURCES = genb_test.c
  genb_test_LDFLAGS = libtsnnls.la 
  TESTS = f77blas_test genb_test $(SHELL_TESTS)
else 
  check_PROGRAMS = f77blas_test transpose_test
  TESTS = f77blas_test $(SHELL_TESTS)
endif

transpose_test_SOURCES = transpose_test.c
transpose_test_LDFLAGS = libtsnnls.la

f77blas_test_SOURCES = f77blas_test.c
