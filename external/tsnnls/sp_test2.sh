#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/A_02.sparse -b $srcdir/b_02.mat -x $srcdir/x_02.mat --tsnnls
else 
./tsnnls_test $srcdir/A_02.sparse $srcdir/b_02.mat $srcdir/x_02.mat --tsnnls
fi
