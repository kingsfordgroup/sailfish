#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/A_01.sparse -b $srcdir/b_01.mat -x $srcdir/x_01.mat --tsnnls --pjv --fallback
else 
./tsnnls_test $srcdir/A_01.sparse $srcdir/b_01.mat $srcdir/x_01.mat --tsnnls
fi
