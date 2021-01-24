#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis

##raw=210105_sindbis
##samp=mAbdpi1
##raw=210107_sindbis
##samp=mAbdpi2
##raw=210112_sindbis
##samp=sinvdpi3

##raw=210114_sindbis
##samp=sinvdpi2

##raw=210116_sindbis
##samp=sinvdpi1

raw=210118_sindbis
samp=mAbdpi3

if [ $1 == org ] ; then
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
    for i in rep1 rep2 rep3 ;
    do
	prefix=${samp}_$i
	rawdir=$datadir/${raw}_$prefix
	mkdir -p $datadir/replicates/$prefix/etc

	cp $rawdir/no_sample/*/*.* $datadir/replicates/$prefix/etc/
    done
fi

if [ $1 == mock ] ; then
    mock=210118_sindbis_mock
    rawdir=$datadir/$mock
    prefix=mock
    
    mkdir -p $datadir/replicates/$prefix

    ##grab fastqs and zip
    mkdir -p $datadir/replicates/$prefix/fqs
    cat $rawdir/no_sample/*/fastq_pass/*fastq > $datadir/replicates/$prefix/fqs/$prefix.fq
    pigz -p 12 $datadir/replicates/$prefix/fqs/$prefix.fq
    
    ##tar the run
    tar czf $datadir/raw/$mock.tar.gz $rawdir
    
    ##upload to aws
    aws s3 cp $datadir/raw/$mock.tar.gz s3://nparchive/sindbis/ --storage-class DEEP_ARCHIVE
fi
