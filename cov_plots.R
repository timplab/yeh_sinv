library(tidyverse)

datadir='/dilithium/Data/Nanopore/sindbis/replicates'
dbxdir='~/Dropbox/timplab_data/sindbis/replicates/cov'

sampinfo=tibble(condition=c(rep('mAb', 9), rep('sinv',9), 'mock'),
                dpi=c(rep(1,3), rep(2,3), rep(3,3), rep(1,3), rep(2,3), rep(3,3), 0),
                rep=c(rep(c(1,2,3), 6), 0))

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
    mutate(max_norm=cov/max(cov))

covplotfile=file.path(dbxdir, 'coverage.pdf')
pdf(covplotfile, w=16, h=9)
for (i in unique(covnorm$cond)) {
    plot=ggplot(covnorm %>% filter(cond==i), aes(x=pos, y=sum_norm, colour=samp)) +
        geom_line() +
        ggtitle(i) +
        xlab('Position') +
        ylab('Normalized Coverage') +
        scale_fill_brewer(palette = "Set2") +
        theme_bw()
    print(plot)
}
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
