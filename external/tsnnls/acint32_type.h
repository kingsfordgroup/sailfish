#ifndef __ACINT32_TYPE_H__
#define __ACINT32_TYPE_H__

/* This file is part of TSNNLS. Here we define a special type which
   should be a 32-bit integer. We use the variables SIZEOF_INT and 
   SIZEOF_LONG_INT, which should be passed in from configure. */

/* We could sidestep this problem completely by using POSIX/C99 types
   such as int32_t and uint32_t, which would probably be the cleanest
   solution, but we don't have an easy way to know where/whether those
   are defined until we upgrade to Autoconf 2.60 or later. */

#ifndef __AC_INT32TYPE_H__
#define __AC_INT32TYPE_H__ 1
   
#include"config.h"
#include "stdint.h"

#if SIZEOF_LONG_INT == 4               

/* This option will be taken in the OSX/Darwin world, where CLAPACK expects
   4-byte (but type long int) integers in calls to the library. */

typedef long int ACINT32_TYPE;
typedef unsigned long int ACUINT32_TYPE;

#elif SIZEOF_INT == 4 

/* We expect this to be the default option in the remaining (sane) world,
   where Lapack's 4-byte integers are just of type int. */

typedef int ACINT32_TYPE;
typedef unsigned int ACUINT32_TYPE;

#else

typedef int32_t ACINT32_TYPE;
typedef uint32_t ACUINT32_TYPE;

//#error Neither int nor long int has size 4. Fix acint32_type.h

#endif
#endif

#endif
