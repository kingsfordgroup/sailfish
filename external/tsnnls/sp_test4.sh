#!/bin/sh
if (test $argtable == 1); then 
./tsnnls_test -A $srcdir/A_04.sparse -b $srcdir/b_04.mat -x $srcdir/x_04.mat --tsnnls
else 
./tsnnls_test $srcdir/A_04.sparse $srcdir/b_04.mat $srcdir/x_04.mat --tsnnls
fi
