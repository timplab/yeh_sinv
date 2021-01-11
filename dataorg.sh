#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis

if [ $1 == org ] ; then

    ##raw=210105_sindbis
    ##samp=mAbdpi1

    raw=210107_sindbis
    samp=mAbdpi2
    
    for i in rep1 rep2 rep3 ;
    do
	prefix=${samp}_$i
	echo $prefix
	rawdir=$datadir/${raw}_$prefix
	mkdir -p $datadir/replicates/$prefix

	##grab fastqs and zip
	mkdir -p $datadir/replicates/$prefix/fqs
	cat $rawdir/no_sample/*/fastq_pass/*fastq > $datadir/replicates/$prefix/fqs/$prefix.fq
	pigz -p 12 $datadir/replicates/$prefix/fqs/$prefix.fq

	##tar the run
	tar czf $datadir/raw/${raw}_$prefix.tar.gz $rawdir

	##upload to aws
	aws s3 cp $datadir/raw/${raw}_$prefix.tar.gz s3://nparchive/sindbis/ --storage-class DEEP_ARCHIVE
    done
fi

	
if [ $1 == grab_misc ] ; then
    raw=210105_sindbis
    samp=mAbdpi1

    for i in rep1 rep2 rep3 ;
    do
	prefix=${samp}_$i
	rawdir=$datadir/${raw}_$prefix
	mkdir -p $datadir/replicates/$prefix/etc

	cp $rawdir/no_sample/*/*.* $datadir/replicates/$prefix/etc/
    done
fi
