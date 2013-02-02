#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/A_05.sparse -b $srcdir/b_05.mat -x $srcdir/x_05.mat --tsnnls
else 
./tsnnls_test $srcdir/A_05.sparse $srcdir/b_05.mat $srcdir/x_05.mat --tsnnls
fi
