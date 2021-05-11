library(tidyverse)
library(RColorBrewer)
library(ggnewscale)
library(cowplot)
library(ShortRead)
library(ggridges)

##datadir='/dilithium/Data/Nanopore/sindbis/replicates'
datadir='/mithril/Data/Nanopore/projects/sindbis/replicates'
dbxdir='~/Dropbox/timplab_data/sindbis/replicates/cov'

sampinfo=tibble(condition=c(rep('mAb', 9), rep('sinv',9), 'mock'),
                dpi=c(rep(1,3), rep(2,3), rep(3,3), rep(1,3), rep(2,3), rep(3,3), 0),
                rep=c(rep(c(1,2,3), 6), 0))
sampinfo$dpi[19]=NA
sampinfo$rep[19]=NA

##annotation
gff=paste0('/dilithium/Data/Nanopore/sindbis/annot_regions.tsv')
regions=read_tsv(gff, col_names=c('prot', 'start', 'end')) %>%
    mutate(colors=brewer.pal(9, 'Set3'))


##get in read count info
aligncountscsv=file.path(dbxdir,'align_counts.csv')
aligncounts=read_csv(aligncountscsv) %>%
    mutate(sinv=total-notsinv)


##read in coverage info
cov_cnames=c('chr','start', 'end','pos','cov')
cov=tibble(chr=as.character(),
           start=as.numeric(),
           end=as.numeric(),
           pos=as.numeric(),
           cov=as.numeric(),
           samp=as.character(),
           cond=as.character(),
           dpi=as.character(),
           rep=as.character())
for (i in 1:dim(sampinfo)[1]) {
    info=sampinfo[i,]
    samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    covfile=file.path(datadir, samp, 'cov', paste0(samp, '.primary.cov'))

    if (file.exists(covfile)) {
        sampcov=read_tsv(covfile, col_names=cov_cnames) %>%
            mutate(samp=samp) %>%
            mutate(cond=paste0(info$condition, 'dpi', as.character(info$dpi))) %>%
            mutate(dpi=info$dpi) %>%
            mutate(rep=info$rep)
        cov=rbind(cov, sampcov)
    }
}
covnorm=cov %>%
    group_by(samp) %>%
    mutate(sum_norm=cov/sum(cov)) %>%
    mutate(max_norm=cov/max(cov)) %>%
    rowwise() %>%
    mutate(num_norm=cov/aligncounts$sinv[aligncounts$samp==samp])

##try with spliced
splicecov=tibble(chr=as.character(),
           start=as.numeric(),
           end=as.numeric(),
           pos=as.numeric(),
           cov=as.numeric(),
           samp=as.character(),
           cond=as.character(),
           dpi=as.character(),
           rep=as.character())
for (i in 1:dim(sampinfo)[1]) {
    info=sampinfo[i,]
    samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    covfile=file.path(datadir, samp, 'cov', paste0(samp, '.splice.cov'))

    if (file.exists(covfile)) {
        sampcov=read_tsv(covfile, col_names=cov_cnames) %>%
            mutate(samp=samp) %>%
            mutate(cond=paste0(info$condition, 'dpi', as.character(info$dpi))) %>%
            mutate(dpi=info$dpi) %>%
            mutate(rep=info$rep)
        splicecov=rbind(splicecov, sampcov)
    }
}
splicenorm=splicecov %>%
    group_by(samp) %>%
    mutate(sum_norm=cov/sum(cov)) %>%
    mutate(max_norm=cov/max(cov)) %>%
    rowwise() %>%
    mutate(num_norm=cov/aligncounts$sinv[aligncounts$samp==samp])


covplotfile=file.path(dbxdir, 'coverage.pdf')
pdf(covplotfile, w=16, h=9)
for (i in unique(covnorm$cond)) {
    plot=ggplot(covnorm %>% filter(cond==i), aes(x=pos, y=sum_norm, colour=samp)) +
        geom_line() +
        ggtitle(i) +
        xlab('Position') +
        ylab('Normalized Coverage') +
        scale_colour_brewer(palette = "Set2") +
        theme_bw()
    print(plot)
}
dev.off()



####averaged across the replicates
meancov=covnorm %>%
    group_by(cond, dpi, pos) %>%
    summarise(avg_numnorm=mean(num_norm)) %>%
    rowwise() %>%
    mutate(status=str_split(cond, 'dpi')[[1]][1]) %>%
    mutate(dpi=as.character(dpi))
excovnorm=covnorm %>%
    filter(samp!='sinvdpi1_rep2')
exmeancov=excovnorm %>%
    group_by(cond, dpi, pos) %>%
    summarise(avg_numnorm=mean(num_norm)) %>%
    rowwise() %>%
    mutate(status=str_split(cond, 'dpi')[[1]][1]) %>%
    mutate(dpi=as.character(dpi))

avgcovplotfile=file.path(dbxdir, 'avg_coverage.pdf')
pdf(avgcovplotfile, w=16, h=9)
ymin=.0001
ymax=.00015
plot=ggplot(meancov, aes(x=pos, y=avg_numnorm, colour=status, linetype=dpi)) +
    geom_line(size=1) +
    ggtitle('Replicate Averages') +
    ylab('Normalized Coverage') +
    scale_colour_brewer(palette = 'Set2') +
    scale_y_log10(limits=c(1e-4,1)) +
    new_scale_color() +
    geom_rect(data=regions, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
    geom_text(data=regions, inherit.aes=FALSE, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot, size=.2)) +
    scale_fill_manual(values=regions$colors, labels=regions$prot) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
print(plot)

explot=ggplot(exmeancov, aes(x=pos, y=avg_numnorm, colour=status, linetype=dpi)) +
    geom_line() +
    ggtitle('Replicate Averages, outlier excluded') +
    ylab('Normalized Coverage') +
    scale_colour_brewer(palette = 'Set2') +
    scale_y_log10() +
    theme_bw()
print(explot)
dev.off()





gencov_cnames=c('chr','pos','cov')
gencov=tibble(chr=as.character(),
           pos=as.numeric(),
           cov=as.numeric(),
           samp=as.character(),
           cond=as.character(),
           dpi=as.character(),
           rep=as.character())

for (i in 1:dim(sampinfo)[1]) {
    info=sampinfo[i,]
    samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    covfile=file.path(datadir, samp, 'cov', paste0(samp, '.primary.genomecov'))

    if (file.exists(covfile)) {
        sampcov=read_tsv(covfile, col_names=gencov_cnames) %>%
            mutate(samp=samp) %>%
            mutate(cond=paste0(info$condition, 'dpi', as.character(info$dpi))) %>%
            mutate(dpi=info$dpi) %>%
            mutate(rep=info$rep)
        gencov=rbind(gencov, sampcov)
    }
}




##plots with genome coverage - means deletions are marked with lower cov. 
gencovnorm=gencov %>%
    group_by(samp) %>%
    mutate(sum_norm=cov/sum(cov)) %>%
    mutate(max_norm=cov/max(cov))

covplotfile=file.path(dbxdir, 'gencoverage.pdf')
pdf(covplotfile, w=16, h=9)
for (i in unique(gencovnorm$cond)) {
    plot=ggplot(gencovnorm %>% filter(cond==i), aes(x=pos, y=sum_norm, colour=samp)) +
        geom_line() +
        ggtitle(i) +
        xlab('Position') +
        ylab('Normalized Coverage') +
        scale_fill_brewer(palette = "Set2") +
        theme_bw()
    print(plot)
}
dev.off()







##make figure
plotcov <- function(covinfo, regions) {
    ##takes in covnorm like object
    ##making colours work: https://stackoverflow.com/questions/58976114/how-do-i-have-multiple-lines-of-the-same-color-with-gg-plot
    mapping=covinfo %>% distinct(cond,samp)
    cols=brewer.pal(6, 'Set2')
    cols=cols[1:n_distinct(mapping$cond)]
    names(cols)=unique(mapping$cond)
    plotcols=cols[mapping$cond]
    names(plotcols)=mapping$samp

    ##barheight=-.07*max(covinfo$sum_norm)
    barheight=-.07*max(covinfo$num_norm)
    regions$ymin=barheight
    regions$ymax=0+.2*barheight
    values=regions$colors
    names(values)=regions$prot
    
    title=as.character(covinfo$dpi[1])
    ##plot=ggplot(covinfo, mapping=aes(x=pos, y=sum_norm, col=samp)) +
    plot=ggplot(covinfo, mapping=aes(x=pos, y=num_norm, col=samp)) +
        geom_line(size=1) +
        ggtitle(paste0('Days post infection: ', title)) +
        xlab('Position') +
        ylab('Normalized Coverage') +
        scale_colour_manual(values = plotcols) +
        new_scale_color() +
        geom_rect(data=regions, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
        geom_text(data=regions, inherit.aes=FALSE, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot, size=.2)) +
        scale_fill_manual(values=regions$colors, labels=regions$prot) +
        theme_bw()
    return(plot)
}

cov1=plotcov(covnorm %>% filter(dpi==1), regions)
cov2=plotcov(covnorm %>% filter(dpi==2), regions)
cov3=plotcov(covnorm %>% filter(dpi==3), regions)


covplotpdf=file.path(dbxdir, 'coverage_fig.pdf')
pdf(covplotpdf, h=18, w=14)
plot_grid(cov1, cov2, cov3, ncol=1, align='v')
dev.off()

covplotpdf_num=file.path(dbxdir, 'coverage_fig_num.pdf')
pdf(covplotpdf_num, h=19, w=16)
plot_grid(cov1, cov2, cov3, ncol=1, align='v')
dev.off()




plotcov_log <- function(covinfo, regions) {
    ##takes in covnorm like object
    ##making colours work: https://stackoverflow.com/questions/58976114/how-do-i-have-multiple-lines-of-the-same-color-with-gg-plot
    mapping=covinfo %>% distinct(cond,samp)
    cols=brewer.pal(6, 'Set2')
    cols=cols[1:n_distinct(mapping$cond)]
    names(cols)=unique(mapping$cond)
    plotcols=cols[mapping$cond]
    names(plotcols)=mapping$samp

    ##barheight=-.07*max(covinfo$sum_norm)
    barheight=-.00025
    regions$ymin=barheight
    regions$ymax=0+.2*barheight
    values=regions$colors
    names(values)=regions$prot

    covinfo=covinfo %>%
        mutate(rep=as.character(rep))
    
    title=as.character(covinfo$dpi[1])
    ##plot=ggplot(covinfo, mapping=aes(x=pos, y=sum_norm, col=samp)) +
    plot=ggplot(covinfo, mapping=aes(x=pos, y=num_norm, col=samp, linetype=rep)) +
        geom_line(size=1) +
        ggtitle(paste0('Days post infection: ', title)) +
        xlab('Position') +
        ylab('Normalized Coverage') +
        scale_colour_manual(values = plotcols) +
        scale_y_log10() +
        theme_bw()
    return(plot)
}

logcov1=plotcov_log(covnorm %>% filter(dpi==1), regions)
logcov2=plotcov_log(covnorm %>% filter(dpi==2), regions)
logcov3=plotcov_log(covnorm %>% filter(dpi==3), regions)

ymin=0
ymax=.1
logrect=ggplot(regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot, alpha=.3)) +
    geom_rect() +
    scale_fill_manual(values=regions$colors, labels=regions$prot) +
    geom_text(data=regions, inherit.aes=FALSE, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot, size=.2)) +
    theme_void()

covplotpdf_num=file.path(dbxdir, 'coverage_fig_lognum.pdf')
pdf(covplotpdf_num, h=19, w=16)
plot_grid(logcov1, logcov2, logcov3, logrect, ncol=1, align='v', rel_heights=c(1,1,1,.06))
dev.off()

spliceplot1=plotcov_log(splicenorm %>% filter(dpi==1), regions)
spliceplot2=plotcov_log(splicenorm %>% filter(dpi==2), regions)
spliceplot3=plotcov_log(splicenorm %>% filter(dpi==3), regions)

spliceplot=file.path(dbxdir, 'coverage_fig_lognum_splice.pdf')
pdf(spliceplot, h=19, w=16)
plot_grid(spliceplot1, spliceplot2, spliceplot3, logrect, ncol=1, align='v', rel_heights=c(1,1,1,.06))
dev.off()






####plot read length distribution

allreadlens=tibble(lengths=as.numeric(),
                   samp=as.character())
for (i in 1:dim(sampinfo)[1]) {
    info=sampinfo[i,]
    if (info$condition=='mock') {
        samp='mock'
    }else{
        samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    }
    print(samp)

    fqfile=file.path(datadir, samp, 'fqs', paste0(samp, '.fq.gz'))
    fqz=gzfile(fqfile, 'rt')
    fq=readLines(fqz)
    seqpos=seq(2, length(fq), 4)
    
    
readlens=tibble(lengths=str_length(fq[seqpos])) %>%
        mutate(samp=samp)
    allreadlens=bind_rows(allreadlens, readlens)
}

interest='sinvdpi1_rep2'
allreadlens=allreadlens %>%
    mutate(status=samp==interest)


pafcols=c('qname', 'qlen', 'qstart', 'qend', 'strand', 'rname', 'rlen', 'rstart', 'rend', 'matches', 'blocklen', 'mapq', 'one', 'two', 'three', 'four', 'five', 'six')
allratlens=tibble(lengths=as.numeric(),
                   samp=as.character())
for (i in 1:dim(sampinfo)[1]) {
    info=sampinfo[i,]
    if (info$condition=='mock') {
        samp='mock'
    }else{
        samp=paste0(info$condition, 'dpi', as.character(info$dpi), '_rep', as.character(info$rep))
    }
    print(samp)

    paffile=file.path(datadir, samp, 'align', paste0(samp, '.rat.paf'))
    lenspaf=read_tsv(paffile, col_names=pafcols) %>%
        filter(mapq>30) %>%
        group_by(qname) %>%
        summarise(lengths=mean(qlen)) %>%
        mutate(samp=samp) %>%
        select(lengths, samp)
    allratlens=bind_rows(allratlens, lenspaf)
}
allratlens=allratlens %>%
    mutate(status=samp==interest)

qcdir='~/Dropbox/timplab_data/sindbis/replicates/qc'
readlenpdf=file.path(qcdir, 'readlens.pdf')
pdf(readlenpdf, w=15, h=8)
plot=ggplot(allreadlens, aes(x=lengths, y=samp, colour=status, fill=status, alpha=.05)) +
    geom_density_ridges() +
    ggtitle('All Read Lengths Log scale') +
    xlab('Read Length') +
    scale_x_log10() +
    theme_bw()
print(plot)
plot=ggplot(allreadlens, aes(x=lengths, y=samp, colour=status, fill=status, alpha=.05)) +
    geom_density_ridges() +
    ggtitle('All Read Lengths Linear scale') +
    xlab('Read Length') +
    xlim(0, 12500) +
    theme_bw()
print(plot)
plot=ggplot(allratlens, aes(x=lengths, y=samp, colour=status, fill=status, alpha=.05)) +
    geom_density_ridges() +
    ggtitle('Rat Read Lengths Log scale') +
    xlab('Read Length') +
    scale_x_log10() +
    theme_bw()
print(plot)
plot=ggplot(allratlens, aes(x=lengths, y=samp, colour=status, fill=status, alpha=.05)) +
    geom_density_ridges() +
    ggtitle('Rat Read Lengths Linear scale') +
    xlim(0, 10000) +
    xlab('Read Length') +
    theme_bw()
print(plot)

dev.off()
