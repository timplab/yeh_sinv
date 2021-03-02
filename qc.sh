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

	
if [ $1 == count_reads ] ; then
    echo samp,numreads >> $dbxdir/qc/readcounts.csv
    linesperread=4
    for samp in $datadir/replicates/* ;
    do
	i=`basename $samp`
	numlines=`zcat $samp/fqs/$i.fq.gz | wc -l `
	echo $i,$((numlines / linesperread)) >> $dbxdir/qc/readcounts.csv
    done
fi
