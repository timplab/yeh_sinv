#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis/replicates
dbxdir=~/Dropbox/timplab_data/sindbis/replicates
ref=/dilithium/Data/Nanopore/sindbis/refs/sindbis_jane.fasta
rat=/dilithium/Data/Nanopore/sindbis/refs/rattus_norvegicus.fa

if [ $1 == align ] ; then
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	mkdir -p $datadir/$i/align
	minimap2 -a -uf -k14 -t 36 $ref $datadir/$i/fqs/$i.fq.gz | \
            samtools view -@ 36 -b | \
            samtools sort -@ 36 -o $datadir/$i/align/$i.sorted.bam
	samtools index $datadir/$i/align/$i.sorted.bam
	
	samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.sorted.bam |
            samtools sort -@ 36 -o $datadir/$i/align/$i.primary.sorted.bam
	samtools index $datadir/$i/align/$i.primary.sorted.bam

	minimap2 -a -uf -k14 -t 36 $rat $datadir/$i/fqs/$i.fq.gz | \
            samtools view -@ 36 -b | \
            samtools sort -@ 36 -o $datadir/$i/align/$i.rat.sorted.bam
	samtools index $datadir/$i/align/$i.rat.sorted.bam
	
	samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.rat.sorted.bam |
            samtools sort -@ 36 -o $datadir/$i/align/$i.rat.primary.sorted.bam
	samtools index $datadir/$i/align/$i.rat.primary.sorted.bam
    done
fi

if [ $1 == cov ] ; then
    mkdir -p $datadir/$i/cov
    bedtools coverage -d -a $ref.bed -b $datadir/$i/align/$i.sorted.bam > $datadir/$i/cov/$i.cov
    bedtools coverage -d -a $ref.bed -b $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.cov
fi

if [ $1 == genomecov ] ; then
    mkdir -p $datadir/$i/cov
    bedtools genomecov -d -ibam $datadir/$i/align/$i.primary.sorted.bam > $datadir/$i/cov/$i.primary.genomecov
fi

if [ $1 == count ] ; then
    echo samp,sinv,rat >> $dbxdir/cov/align_counts.csv
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	sinv=`samtools view -c -F 260 $samp/align/$i.primary.sorted.bam`
	rat=`samtools view -c -F 260 $samp/align/$i.rat.primary.sorted.bam`
	echo $i,$sinv,$rat >> $dbxdir/cov/align_counts.csv
    done
fi
