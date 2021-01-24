#!/bin/bash

srcdir=~/Code/yfan_nanopore/sindbis/replicates

if [ $1 == set1 ] ; then
    for i in mAbdpi1_rep1 mAbdpi1_rep2 mAbdpi1_rep3 ;
    do
	bash $srcdir/coverage.sh align $i
	bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi

if [ $1 == set2 ] ; then
    for i in mAbdpi2_rep1 mAbdpi2_rep2 mAbdpi2_rep3 ;
    do
	bash $srcdir/coverage.sh align $i
	bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi

if [ $1 == set3 ] ; then
    for i in sinvdpi3_rep1 sinvdpi3_rep2 sinvdpi3_rep3 ;
    do
	bash $srcdir/coverage.sh align $i
	bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi

if [ $1 == set4 ] ; then
    for i in sinvdpi2_rep1 sinvdpi2_rep2 sinvdpi2_rep3 ;
    do
	bash $srcdir/coverage.sh align $i
	bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi

if [ $1 == set5 ] ; then
    for i in sinvdpi1_rep1 sinvdpi1_rep2 sinvdpi1_rep3 ;
    do
	bash $srcdir/coverage.sh align $i
	bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi

if [ $1 == set6 ] ; then
    for i in mAbdpi3_rep1 mAbdpi3_rep2 mAbdpi3_rep3 ;
    do
	bash $srcdir/coverage.sh align $i
	bash $srcdir/coverage.sh cov $i
	bash $srcdir/coverage.sh genomecov $i
    done
fi
	
