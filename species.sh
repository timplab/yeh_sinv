#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/sindbis/replicates
dbxdir=~/Dropbox/timplab_data/sindbis/replicates
#ref=/mithril/Data/Nanopore/sindbis/refs/sindbis_jane.fasta
ref=/mithril/Data/Nanopore/projects/sindbis/refs/sindbis_jane_species.fasta
rat=/mitrhil/Data/Nanopore/sindbis/refs/rattus_norvegicus.fa

if [ $1 == align_species ] ; then
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	mkdir -p $datadir/$i/align
	minimap2 -uf -k14 -t 36 $ref $datadir/$i/fqs/$i.fq.gz > \
                 $datadir/$i/align/$i.species.paf
		
    done
fi
	
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
    echo samp,notsinv,notrat,total >> $dbxdir/cov/align_counts.csv
    for samp in $datadir/* ;
    do
	i=`basename $samp`

	##count reads that did not align
	notsinv=`samtools view -c -f 4 $samp/align/$i.sorted.bam`
	notrat=`samtools view -c -f 4 $samp/align/$i.rat.sorted.bam`

	
	fqlines=`zcat $samp/fqs/$i.fq.gz | wc -l`
	fqdiv=4
	total=`echo $((fqlines / fqdiv))`
	
	echo $i,$notsinv,$notrat,$total >> $dbxdir/cov/align_counts.csv
    done
fi


if [ $1 == splice_align ] ; then
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	minimap2 -ax splice -uf -k14 --splice-flank=no -t 36 $ref $datadir/$i/fqs/$i.fq.gz | 
	    samtools view -@ 36 -b | 
	    samtools sort -@ 36 -o $datadir/$i/align/$i.splice.sorted.bam
	samtools index $datadir/$i/align/$i.splice.sorted.bam
    done
fi
     
if [ $1 == splice_filt ] ; then
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.splice.sorted.bam |
	    samtools sort -@ 36 -o $datadir/$i/align/$i.splice.primary.sorted.bam

    done
fi

if [ $1 == splice_cov ] ; then
    for samp in $datadir/* ;
    do
	(i=`basename $samp`
	
        bedtools coverage -d -split \
		 -a $ref.bed \
		 -b $datadir/$i/align/$i.splice.sorted.bam > $datadir/$i/cov/$i.splice.cov
	bedtools coverage -d -split \
		 -a $ref.bed \
		 -b $datadir/$i/align/$i.splice.primary.sorted.bam > $datadir/$i/cov/$i.splice.primary.cov ) &
    done
fi


if [ $1 == stringtie ] ; then
    for samp in $datadir/* ;
    do
	i=`basename $samp`
	mkdir -p $datadir/$i/stringtie

	stringtie \
	    $datadir/$i/align/$i.splice.sorted.bam \
	    -o $datadir/$i/stringtie/$i.splice.gff \
	    -L \
	    -p 36
    done
fi

	
	
       
if [ $1 == stringtie_vanilla ] ; then
    for samp in $datadir/* ;
    do
        i=`basename $samp`
        mkdir -p $datadir/$i/stringtie
	stringtie \
            $datadir/$i/align/$i.sorted.bam \
            -o $datadir/$i/stringtie/$i.gff \
            -L \
            -p 36
    done
fi



if [ $1 == pafalign ] ; then
    for samp in $datadir/* ;
    do
        i=`basename $samp`
        minimap2 -uf -k14 -t 36 $rat $datadir/$i/fqs/$i.fq.gz > \
                 $datadir/$i/align/$i.rat.paf
    done
fi
