#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/A_06.sparse -b $srcdir/b_06.mat -x $srcdir/x_06.mat --tsnnls
else 
./tsnnls_test $srcdir/A_06.sparse $srcdir/b_06.mat $srcdir/x_06.mat --tsnnls
fi
