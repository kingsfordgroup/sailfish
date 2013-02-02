#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/rrtestA.sparse -b $srcdir/rrtestb.mat -x $srcdir/rrtestx.mat --tsnnls
else 
./tsnnls_test $srcdir/rrtestA.sparse $srcdir/rrtestb.mat $srcdir/rrtestx.mat --tsnnls
fi
