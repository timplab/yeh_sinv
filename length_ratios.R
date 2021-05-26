library(tidyverse)
library(GenomicRanges)
library(cowplot)

datadir='/mithril/Data/Nanopore/projects/sindbis/replicates'
dbxdir='~/Dropbox/timplab_data/sindbis/replicates/cov'
pafcols=c('qname', 'qlen', 'qstart', 'qend', 'starnd', 'rname', 'rlen', 'rstart', 'rend', 'matches', 'alen', 'mapq', 'type', 'cm', 's1', 's2', 'dv', 'rl')

sampinfo=tibble(condition=c(rep('mAb', 9), rep('sinv',9), 'mock'),
                dpi=c(rep(1,3), rep(2,3), rep(3,3), rep(1,3), rep(2,3), rep(3,3), 0),
                rep=c(rep(c(1,2,3), 6), 0))
sampinfo$dpi[19]=NA
sampinfo$rep[19]=NA


####find drop points using the cov data
##going through partially by hand to figure out best cutoffs to classify dvg vs subgenomic, etc
covfile=file.path(datadir, 'sinvdpi1_rep1', 'cov', 'sinvdpi1_rep1.primary.cov')
cov_cols=c('chr','start', 'end','pos','cov')
cov=read_tsv(covfile, col_names=cov_cols)
deltacov=diff(cov$cov)

regions=GRanges(seqnames=c('dvg', 'subdvg', 'subgenomic'),
                ranges=IRanges(start=c(10, 10, 7609),
                               end=c(1140, 1750, 11700)))


ratios=tibble(name=as.character(),
              subgen=as.numeric(),
              dvg=as.numeric(),
              subdvg=as.numeric(),
              condition=as.character(),
              dpi=as.numeric())

for (i in 1:(dim(sampinfo)[1]-1)) {
    info=sampinfo[i,]
    samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    paffile=file.path(datadir, samp, 'align', paste0(samp, '.paf'))

    lenspaf=read_tsv(paffile, col_names=pafcols) %>%
        filter(mapq>20) %>%
        select(qname, qlen, rname, rlen, rstart, rend, mapq, type)
    full=lenspaf %>%
        filter(rstart<7000) %>%
        filter(rend>3000)
    #defining drops are 7608 to 11694
    subgen=lenspaf %>%
        filter(rstart>7610 & rstart<11690)
    ##dvg drop at 1751
    dvg=lenspaf %>%
        filter(rend>1749 & rend<1755)
    subdvg=lenspaf %>%
        ##filter(qlen>100) %>%
        ##filter(rstart<50) %>%
        filter(rend>50 & rend<1745)

    subgen_ratio=dim(subgen)[1]/dim(full)[1]
    dvg_ratio=dim(dvg)[1]/dim(full)[1]
    subdvg_ratio=dim(subdvg)[1]/dim(full)[1]

    ratioinfo=tibble(name=samp,
                     subgen=subgen_ratio,
                     dvg=dvg_ratio,
                     subdvg=subdvg_ratio,
                     condition=info$condition,
                     dpi=info$dpi)
    ratios=bind_rows(ratios, ratioinfo)
}

ratios$dpi=as.character(ratios$dpi)

speciespdf=file.path(dbxdir, 'rna_species_ratios.pdf')
pdf(speciespdf, h=5, w=13)
dvgplot=ggplot(ratios, aes(x=dpi, y=dvg, colour=condition, fill=condition, alpha=.2)) +
    geom_point(size=6) +
    ggtitle('DVG') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
subdvgplot=ggplot(ratios, aes(x=dpi, y=subdvg, colour=condition, fill=condition, alpha=.2)) +
    geom_point(size=6) +
    ggtitle('Sub DVG') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
subgenplot=ggplot(ratios, aes(x=dpi, y=subgen, colour=condition, fill=condition, alpha=.2)) +
    geom_point(size=6) +
    ggtitle('Subgenomic') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
plot_grid(dvgplot, subdvgplot, subgenplot, ncol=3)

dvgplot=ggplot(ratios, aes(x=dpi, y=dvg, colour=condition, fill=condition, alpha=.2)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('DVG') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
subdvgplot=ggplot(ratios, aes(x=dpi, y=subdvg, colour=condition, fill=condition, alpha=.2)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('Sub DVG') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
subgenplot=ggplot(ratios, aes(x=dpi, y=subgen, colour=condition, fill=condition, alpha=.2)) +
    geom_jitter(size=6, width=.05) +
    ggtitle('Subgenomic') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
plot_grid(dvgplot, subdvgplot, subgenplot, ncol=3)
dev.off()

speciespdf=file.path(dbxdir, 'rna_species_ratios_nooutlier.pdf')
pdf(speciespdf, h=5, w=12)
ggplot(ratios, aes(x=dpi, y=dvg, colour=condition, fill=condition, alpha=.2)) +
    geom_point(size=6) +
    ggtitle('DVG') +
    ylim(0,1) +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
ggplot(ratios, aes(x=dpi, y=subdvg, colour=condition, fill=condition, alpha=.2)) +
    geom_point(size=6) +
    ggtitle('Sub DVG') +
    ylim(0,4) +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
ggplot(ratios, aes(x=dpi, y=subgen, colour=condition, fill=condition, alpha=.2)) +
    geom_point(size=6) +
    ggtitle('Subgenomic') +
    ylim(0,250) +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
dev.off()






####ratios based on total reads as suggested by jane
ratios=tibble(name=as.character(),
              full=as.numeric(),
              subgen=as.numeric(),
              dvg=as.numeric(),
              subdvg=as.numeric(),
              condition=as.character(),
              dpi=as.numeric())
rationums=tibble(name=as.character(),
              full=as.numeric(),
              subgen=as.numeric(),
              dvg=as.numeric(),
              subdvg=as.numeric(),
              condition=as.character(),
              dpi=as.numeric(),
              total=as.numeric())

for (i in 1:(dim(sampinfo)[1]-1)) {
    info=sampinfo[i,]
    samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    paffile=file.path(datadir, samp, 'align', paste0(samp, '.paf'))
    
    lenspaf=read_tsv(paffile, col_names=pafcols) %>%
        filter(mapq>20) %>%
        select(qname, qlen, rname, rlen, rstart, rend, mapq, type)
    totlen=dim(lenspaf)[1]
    
    full=lenspaf %>%
        filter(rstart<7000) %>%
        filter(rend>3000)
    totfull=dim(full)[1]
    ##defining drops are 7608 to 11694
    subgen=lenspaf %>%
        filter(rstart>7610 & rstart<11690)
    totsubgen=dim(subgen)[1]
    ##dvg drop at 1751
    dvg=lenspaf %>%
        filter(rend>1749 & rend<1755)
    totdvg=dim(dvg)[1]
    subdvg=lenspaf %>%
        ##filter(qlen>100) %>%
        ##filter(rstart<50) %>%
        filter(rend>50 & rend<1745)
    totsubdvg=dim(subdvg)[1]
    
    full_ratio=totfull/totlen
    subgen_ratio=totsubgen/totlen
    dvg_ratio=totdvg/totlen
    subdvg_ratio=totsubdvg/totlen
    
    ratioinfo=tibble(name=samp,
                     full=full_ratio,
                     subgen=subgen_ratio,
                     dvg=dvg_ratio,
                     subdvg=subdvg_ratio,
                     condition=info$condition,
                     dpi=info$dpi)
    ratios=bind_rows(ratios, ratioinfo)
    numinfo=tibble(name=samp,
                   full=totfull,
                   subgen=totsubgen,
                   dvg=totdvg,
                   subdvg=totsubdvg,
                   condition=info$condition,
                   dpi=info$dpi,
                   total=totlen)
    rationums=bind_rows(rationums, numinfo)
}

rationums=rationums %>%
    select(-condition, -dpi)
spcountsfile=file.path(dbxdir, 'ratio_counts.csv')
write_csv(rationums, spcountsfile)

ratios$dpi=as.character(ratios$dpi)
plot_ratios <- function(df, name) {
    plot=ggplot(df, aes(x=dpi, y=ratio, colour=condition, fill=condition)) +
        geom_point(size=6, position=position_dodge(width=.3), alpha=.3) +
        ggtitle(name) +
        stat_summary(fun.y='mean', geom='crossbar', 
                     mapping=aes(ymin=..y.., ymax=..y..), width=.33,
                     show.legend = FALSE) +
        scale_fill_brewer(palette='Set2') +
        scale_colour_brewer(palette='Set2') +
        theme_bw() +
        theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
    return(plot)
}
plot_ratios_aligned <- function(df, name) {
    plot=ggplot(df, aes(x=dpi, y=ratio, colour=condition, fill=condition)) +
        geom_point(size=6, alpha=.3) +
        ggtitle(name) +
        stat_summary(fun.y='mean', geom='crossbar', 
                     mapping=aes(ymin=..y.., ymax=..y..), width=.33,
                     show.legend = FALSE) +
        scale_fill_brewer(palette='Set2') +
        scale_colour_brewer(palette='Set2') +
        theme_bw() +
        theme(panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
    return(plot)
}

full=ratios %>%
    rename('full'='ratio')
subgen=ratios %>%
    rename('subgen'='ratio')
dvg=ratios %>%
    rename('dvg'='ratio')
subdvg=ratios %>%
    rename('subdvg'='ratio')

speciesyieldpdf=file.path(dbxdir, 'rna_species_ratios_numreads.pdf')
pdf(speciesyieldpdf, h=5, w=17)
fullplot=plot_ratios(full, 'Genomic')
dvgplot=plot_ratios(dvg, 'DVG')
subdvgplot=plot_ratios(subdvg, 'sub DVG')
subgenplot=plot_ratios(subgen, 'Subgenomic')
plot_grid(fullplot, dvgplot, subdvgplot, subgenplot, ncol=4)
dev.off()

speciesyieldalignedpdf=file.path(dbxdir, 'rna_species_ratios_numreads_aligned.pdf')
pdf(speciesyieldalignedpdf, h=5, w=17)
fullplot=plot_ratios_aligned(full, 'Genomic')
dvgplot=plot_ratios_aligned(dvg, 'DVG')
subdvgplot=plot_ratios_aligned(subdvg, 'sub DVG')
subgenplot=plot_ratios_aligned(subgen, 'Subgenomic')
plot_grid(fullplot, dvgplot, subdvgplot, subgenplot, ncol=4)
dev.off()


ratioscsv=file.path(dbxdir, 'ratios.csv')
write_csv(ratios, ratioscsv)





