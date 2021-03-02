library(tidyverse)

dbxdir='~/Dropbox/timplab_data/sindbis/replicates'
aligncountscsv=file.path(dbxdir, 'cov','align_counts.csv')

aligncounts=read_csv(aligncountscsv) %>%
    mutate(persinv=sinv/(sinv+rat)) %>%
    mutate(perrat=rat/(sinv+rat)) %>%
    select(-sinv, -rat) %>%
    gather('species', 'value', -samp)

aligncountspdf=file.path(dbxdir, 'cov', 'align_counts.pdf')
pdf(aligncountspdf, h=9, w=16)
ggplot(aligncounts, aes(x=samp, y=value, colour=species, fill=species, alpha=.7)) +
    geom_bar(position='stack', stat='identity') +
    ggtitle('Read Ratios') +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
dev.off()
   
