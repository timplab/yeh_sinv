#!/bin/bash

srcdir=~/Code/yfan_nanopore/sindbis/replicates

if [ $1 == set1 ] ; then

    for i in mAbdpi1_rep1 mAbdpi1_rep2 mAbdpi1_rep3 ;
    do
	#bash $srcdir/coverage.sh align $i
	##bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi
