#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis

if [ $1 == untar ] ; then
    mkdir -p ~/data/sindbis
    mkdir -p ~/data/sindbis/raw
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	( prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	mkdir -p ~/data/sindbis/raw/$prefix
	tar -xzf $i -C ~/data/sindbis/raw/$prefix ) & 
    done
fi

if [ $1 == multi ] ; then
    mkdir -p ~/data/sindbis/multiraw
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	mkdir -p ~/data/sindbis/multiraw/$prefix
	single_to_multi_fast5 \
	    -i ~/data/sindbis/raw/$prefix \
	    -s ~/data/sindbis/multiraw/$prefix \
	    -f $prefix \
	    -n 4000 \
	    -t 36 \
	    --recursive
    done
fi

if [ $1 == call ] ; then
    mkdir -p ~/data/sindbis/called
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	mkdir -p ~/data/sindbis/called/$prefix
	guppy_basecaller \
	    --recursive \
	    --compress_fastq \
	    --kit SQK-RNA001 \
	    --flowcell FLO-MIN106 \
	    -i ~/data/sindbis/multiraw/$prefix \
	    -s ~/data/sindbis/called/$prefix \
	    -x 'cuda:0'
    done
fi

if [ $1 == gather ] ; then
    mkdir -p ~/data/sindbis/fqs
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	( prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	  mkdir -p ~/data/sindbis/fqs/$prefix
	  cat ~/data/sindbis/called/$prefix/*fastq.gz > ~/data/sindbis/fqs/$prefix/${prefix}_multiraw.fq.gz ) &
    done
fi

if [ $1 == call_singlefast5 ] ; then
    ##saw some weird error messages when converting to multifast5 checking single fast5 too
    mkdir -p ~/data/sindbis/called_single
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	mkdir -p ~/data/sindbis/called_single/$prefix
	guppy_basecaller \
	    --recursive \
	    --compress_fastq \
	    --kit SQK-RNA001 \
	    --flowcell FLO-MIN106 \
	    -i ~/data/sindbis/raw/$prefix \
	    -s ~/data/sindbis/called_single/$prefix \
	    -x 'cuda:0'
    done
fi

if [ $1 == gather_singlefast5 ] ; then
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	( prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	  mkdir -p ~/data/sindbis/fqs/$prefix
	  cat ~/data/sindbis/called_single/$prefix/*fastq.gz > ~/data/sindbis/fqs/$prefix/${prefix}.fq.gz ) &
    done
fi

if [ $1 == copy ] ; then
    for i in $datadir/raw/18*TE*.tar.gz ;
    do
	( prefix=`basename $i .tar.gz | cut -d _ -f 2-`
	  mkdir -p $datadir/singlet/$prefix/fqs
	  cp ~/data/sindbis/fqs/$prefix/$prefix.fq.gz $datadir/singlet/$prefix/fqs/$prefix.fq.gz ) &
    done
fi

ref=/dilithium/Data/Nanopore/sindbis/refs/sindbis_jane.fasta
rat=/dilithium/Data/Nanopore/sindbis/refs/rattus_norvegicus.fa

if [ $1 == align ] ; then
    for i in $datadir/singlet/* ;
    do
	prefix=`basename $i`
	mkdir -p $i/align
	
	minimap2 -a -uf -k14 -t 36 $ref $i/fqs/$prefix.fq.gz | \
            samtools view -@ 36 -b | \
            samtools sort -@ 36 -o $i/align/$prefix.sorted.bam
	samtools index $i/align/$prefix.sorted.bam
	
	samtools view -@ 36 -b -F 0x100 $i/align/$prefix.sorted.bam |
            samtools sort -@ 36 -o $i/align/$prefix.primary.sorted.bam
	samtools index $i/align/$prefix.primary.sorted.bam

	minimap2 -a -uf -k14 -t 36 $rat $i/fqs/$prefix.fq.gz | \
            samtools view -@ 36 -b | \
            samtools sort -@ 36 -o $i/align/$prefix.rat.sorted.bam
	samtools index $i/align/$prefix.rat.sorted.bam
	
	samtools view -@ 36 -b -F 0x100 $i/align/$prefix.rat.sorted.bam |
            samtools sort -@ 36 -o $i/align/$prefix.rat.primary.sorted.bam
	samtools index $i/align/$prefix.rat.primary.sorted.bam
    done
fi

if [ $1 == cov ] ; then
    for i in $datadir/singlet/* ;
    do
	prefix=`basename $i`
	mkdir -p $i/cov
	bedtools coverage -d -a $ref.bed -b $i/align/$prefix.sorted.bam > $i/cov/$prefix.cov
	bedtools coverage -d -a $ref.bed -b $i/align/$prefix.primary.sorted.bam > $i/cov/$prefix.primary.cov
    done
fi

if [ $1 == genomecov ] ; then
    for samp in $datadir/singlet/* ;
    do
	i=`basename $samp`
	mkdir -p $i/cov
	bedtools genomecov -d -ibam $i/align/$prefix.primary.sorted.bam > $i/cov/$prefix.primary.genomecov
    done
fi
