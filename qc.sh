#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis
dbxdir=~/Dropbox/timplab_data/sindbis/replicates

if [ $1 == set1_runsum ] ; then

    for i in 1 2 3 ;
    do
	prefix=mAbdpi1_rep$i
	sumfile=$datadir/210105_sindbis_mAbdpi1_rep$i/no_sample/*/sequencing_summary*.txt
	Rscript ~/Code/utils/qc/run_summary.R \
		-i $sumfile \
		-o $dbxdir/$prefix.pdf \
		-p $prefix
    done
fi

	
