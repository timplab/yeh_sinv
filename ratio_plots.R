library(tidyverse)
library(cowplot)

dbxdir='~/Dropbox/timplab_data/sindbis/replicates'

aligncountscsv=file.path(dbxdir, 'cov','align_counts.csv')
aligninfo=read_csv(aligncountscsv)

aligncounts=read_csv(aligncountscsv) %>%
    mutate(sinv=total-notsinv) %>%
    mutate(rat=total-notrat) %>%
    mutate(persinv=sinv/(sinv+rat)) %>%
    mutate(perrat=rat/(sinv+rat)) %>%
    mutate(pertotsinv=sinv/total) %>%
    mutate(pertotsrat=rat/total) #%>%
    select(-sinv, -rat, -notrat, -notsinv, -total) %>%
    gather('species', 'value', -samp) %>%
    rowwise() %>%
    mutate(group=str_split(samp, '_')[[1]][1]) %>%
    mutate(day=str_split(group, 'dpi')[[1]][2])

aligncountspdf=file.path(dbxdir, 'cov', 'align_counts.pdf')
pdf(aligncountspdf, h=9, w=16)
ggplot(aligncounts, aes(x=samp, y=value, colour=species, fill=species, alpha=.7)) +
    geom_bar(position='stack', stat='identity') +
    ggtitle('Read Ratios') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
dev.off()

pointcounts=aligncounts %>%
    filter(species=='persinv')
aligncountpointspdf=file.path(dbxdir, 'cov', 'align_counts_points.pdf')
pdf(aligncountpointspdf, h=9, w=11)
ggplot(pointcounts, aes(x=group, y=value, colour=species, fill=species, alpha=.3)) +
    geom_point(size=5) +
    ggtitle('Read Ratios') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
dev.off()


plot_dayratio <- function(info) {
    ##ratio point plots
    ##take subset of aligncounts from above
    title=paste0('Days post infection: ', info$day[1])
    plot=ggplot(info, aes(x=group, y=value, colour=group, fill=group, alpha=.7)) +
        geom_point(size=6) +
        ggtitle(title) +
        ylim(0,1) +
        scale_fill_brewer(palette='Set2') +
        scale_colour_brewer(palette='Set2') +
        theme_bw()
    return(plot)
}


sinvcounts=aligncounts %>%
    filter(species=='persinv')
day1plot=plot_dayratio(sinvcounts %>% filter(day==1))
day2plot=plot_dayratio(sinvcounts %>% filter(day==2))
day3plot=plot_dayratio(sinvcounts %>% filter(day==3))
mockplot=plot_dayratio(sinvcounts %>% filter(group=='mock'))
ratiospdf=file.path(dbxdir, 'cov', 'aligncounts_fig.pdf')
pdf(ratiospdf, h=5, w=18)
plot_grid(mockplot, day1plot, day2plot, day3plot, ncol=4, align='h')
dev.off()

###check singlet
alignsingletcsv=file.path(dbxdir, 'cov', 'align_counts_singlet.csv')
alignsinglet=read_csv(alignsingletcsv) %>%
    mutate(persinv=sinv/(sinv+rat)) %>%
    mutate(perrat=rat/(sinv+rat)) %>%
    select(-sinv, -rat) %>%
    gather('species', 'value', -samp)

alignsingletpdf=file.path(dbxdir, 'cov', 'align_counts_singlet.pdf')
pdf(alignsingletpdf, h=9, w=16)
ggplot(alignsinglet, aes(x=samp, y=value, colour=species, fill=species, alpha=.7)) +
    geom_bar(position='stack', stat='identity') +
    ggtitle('Read Ratios') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
dev.off()
