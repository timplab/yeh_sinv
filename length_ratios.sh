#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/sindbis/replicates
ref=/mithril/Data/Nanopore/projects/sindbis/refs/sindbis_jane.fasta

if [ $1 == getpaf ] ; then
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	minimap2 -uf -k14 -t 36 $ref $datadir/$i/fqs/$i.fq.gz > \
		 $datadir/$i/align/$i.paf
    done
fi

