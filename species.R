library(tidyverse)
library(cowplot)


##same species ratio analysis as length_ratios.R but based on alignment to transcriptome
datadir='/mithril/Data/Nanopore/projects/sindbis/replicates'
dbxdir='~/Dropbox/timplab_data/sindbis/replicates/cov'
pafcols=c('qname', 'qlen', 'qstart', 'qend', 'starnd', 'rname', 'rlen', 'rstart', 'rend', 'matches', 'alen', 'mapq', 'type', 'cm', 's1', 's2', 'dv', 'rl')

sampinfo=tibble(condition=c(rep('mAb', 9), rep('sinv',9), 'mock'),
                dpi=c(rep(1,3), rep(2,3), rep(3,3), rep(1,3), rep(2,3), rep(3,3), 0),
                rep=c(rep(c(1,2,3), 6), 0))
sampinfo$dpi[19]=NA
sampinfo$rep[19]=NA

readcounts=tibble(sample=as.character(),
                  genomic=as.numeric(),
                  subgenomic=as.numeric(),
                  dvg=as.numeric(),
                  subdvg=as.numeric(),
                  total=as.numeric(),
                  condition=as.character(),
                  dpi=as.numeric())

for (i in 1:(dim(sampinfo)[1]-1)) {
    info=sampinfo[i,]
    samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    paffile=file.path(datadir, samp, 'align', paste0(samp, '.species.paf'))

    lenspaf=read_tsv(paffile, col_names=pafcols) %>%
        select(qname, qlen, rname, rlen, rstart, rend, mapq, type)

    genomic=lenspaf %>%
        filter(rname=='genomic') %>%
        filter(qlen>11000)
    subgenomic=lenspaf %>%
        filter(rname=='subgenomic') %>%
        filter(qlen>4000) %>%
        filter(qlen<4100)
    dvg=lenspaf %>%
        filter(rname=='dvg') %>%
        filter(qlen>1680) %>%
        filter(qlen<1700)
    subdvg=lenspaf %>%
        filter(rname=='subdvg') %>%
        filter(qlen>850) %>%
        filter(qlen<950)

    counts=tibble(sample=samp, genomic=dim(genomic)[1], subgenomic=dim(subgenomic)[1], dvg=dim(dvg)[1], subdvg=dim(subdvg)[1], total=dim(lenspaf)[1], condition=info$condition, dpi=info$dpi)
    readcounts=bind_rows(readcounts, counts)
}

readcounts$dpi=as.character(readcounts$dpi)

ratios=readcounts %>%
    mutate(pergen=genomic/total) %>%
    mutate(persubgen=subgenomic/total) %>%
    mutate(perdvg=dvg/total) %>%
    mutate(persubdvg=subdvg/total) %>%
    select(-genomic, -subgenomic, -dvg, -subdvg)

speciesyieldpdf=file.path(dbxdir, 'rna_species_ratios_numreads_alignment.pdf')
pdf(speciesyieldpdf, h=5, w=17)
fullplot=ggplot(ratios, aes(x=dpi, y=pergen, colour=condition, fill=condition, alpha=.1)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('Genomic RNA') +
    stat_summary(fun.y='mean', geom='crossbar',
                 mapping=aes(ymin=..y.., ymax=..y..), width=.5,
                 show.legend = FALSE) +
            scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
subplot=ggplot(ratios, aes(x=dpi, y=persubgen, colour=condition, fill=condition, alpha=.1)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('Subgenomic RNA') +
    stat_summary(fun.y='mean', geom='crossbar',
                 mapping=aes(ymin=..y.., ymax=..y..), width=.5,
                 show.legend = FALSE) +
            scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
dvgplot=ggplot(ratios, aes(x=dpi, y=perdvg, colour=condition, fill=condition, alpha=.1)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('dvg RNA') +
    stat_summary(fun.y='mean', geom='crossbar',
                 mapping=aes(ymin=..y.., ymax=..y..), width=.5,
                 show.legend = FALSE) +
            scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
subdvgplot=ggplot(ratios, aes(x=dpi, y=persubdvg, colour=condition, fill=condition, alpha=.1)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('subdvg RNA') +
    stat_summary(fun.y='mean', geom='crossbar',
                 mapping=aes(ymin=..y.., ymax=..y..), width=.5,
                 show.legend = FALSE) +
            scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
plot_grid(fullplot, subplot, dvgplot, subdvgplot, ncol=4)
dev.off()



####genomic vs subgenomic ratios
highcountscsv=file.path(dbxdir, 'ratio_counts.csv')
highcounts=read_csv(highcountscsv)
high_est=highcounts %>%
    mutate(ratio=subgen/full) %>%
    filter(dpi==2)

low_est=readcounts %>%
    mutate(ratio=subgenomic/genomic)
bardata=low_est %>%
    group_by(condition, dpi) %>%
    summarise(meanratio=mean(ratio), sdratio=sd(ratio))
    
bars2dpipdf=file.path(dbxdir, 'ratio_counts.pdf')
pdf(bars2dpipdf)
ggplot(bardata, aes(x=dpi, y=meanratio, colour=condition, fill=condition, alpha=.5)) +
    geom_bar(stat='identity', position=position_dodge()) +
    geom_errorbar(aes(ymin=meanratio, ymax=meanratio+sdratio), position=position_dodge(.9), width=.2) +
    ggtitle('subgenomic/genomic') +
    theme_bw()
dev.off()

low_test=low_est %>%
    filter(dpi=='2')
a=t.test(low_test$ratio[1:3], low_test$ratio[4:6])
