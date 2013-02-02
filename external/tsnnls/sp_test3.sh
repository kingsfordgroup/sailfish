#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/A_03.sparse -b $srcdir/b_03.mat -x $srcdir/x_03.mat --tsnnls
else 
./tsnnls_test $srcdir/A_03.sparse $srcdir/b_03.mat $srcdir/x_03.mat --tsnnls
fi
