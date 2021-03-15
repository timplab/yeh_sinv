library(tidyverse)

dbxdir='~/Dropbox/timplab_data/sindbis/replicates'

aligncountscsv=file.path(dbxdir, 'cov','align_counts.csv')
aligncounts=read_csv(aligncountscsv) %>%
    mutate(persinv=sinv/(sinv+rat)) %>%
    mutate(perrat=rat/(sinv+rat)) %>%
    select(-sinv, -rat) %>%
    gather('species', 'value', -samp) %>%
    rowwise() %>%
    mutate(group=str_split(samp, '_')[[1]][1])

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
ggplot(pointcounts, aes(x=group, y=value, colour=species, fill=species, alpha=.7)) +
    geom_point(size=5) +
    ggtitle('Read Ratios') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
dev.off()
   



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
