#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script cuts the appropriate columns from the expression file
and aggregates expression esitmates to the gene level

OPTIONS:
   -h      Show this message
   -c      List of columns to extract
   -n      Normalization method to perform when aggregating
   -i      Input file
EOF
}

lines=
method=
init=
while getopts “hc:n:i:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         c)
             lines=$OPTARG
             ;;
         n)
             method=$OPTARG
             ;;
         i)
             init=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


echo "cut -f ${lines} ${init} | tail -n +2 > ${init}.cut"
cut -f ${lines} ${init} | tail -n +2 > ${init}.cut

echo "python AggregateToGeneLevel.py --gtf=f_NM_ONLY_hg18.refseq.gtf --out=${init}.cut.genes --texp=${init}.cut --method=${method}"
python AggregateToGeneLevel.py --gtf=f_NM_ONLY_hg18.refseq.gtf --out=${init}.cut.genes --texp=${init}.cut --method=${method}