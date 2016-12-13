#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(readr)

otu_table = read_tsv('../../data/rdp_g.melt') %>%
  # compute relative abundances
  group_by(sample) %>%
  mutate(sample_counts=sum(counts), ra=counts/sample_counts) %>%
  ungroup %>%
  # remove samples with few counts
  filter(sample_counts >= 1000)

plot_timeseries = function(this.otu) {
  # create two separate data sets for plotting
  all_data = otu_table %>% filter(otu == this.otu, day <= 17)
  hsd_data = all_data %>% filter(diet == 'H')
  nd_data = all_data %>% filter(diet == 'N')

  all_data %>%
    ggplot(aes(x=factor(day), y=ra, color=diet)) +
      stat_boxplot(data=hsd_data, outlier.shape=NA, color='orange', fill=NA) +
      stat_boxplot(data=nd_data, outlier.shape=NA, color='black', fill=NA) +
      geom_point(data=hsd_data, color='orange', position=position_jitter(h=0, w=0.1)) +
      geom_point(data=nd_data, color='black', position=position_jitter(h=0, w=0.1)) +
      xlab('day') +
      ylab('relative abundance') +
      theme_minimal()
}

plot_timeseries_ymax = function(this.otu, ymax) {
  plot_timeseries(this.otu) + ylim(0.0, ymax)
}

plot_timeseries_day14_ymax = function(this.otu) {
  ymax = otu_table %>% filter(otu == this.otu, day == 14) %$% max(ra)
  plot_timeseries_ymax(this.otu, ymax)
}

otus = c('Root;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus',
         'Root;Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;Parasutterella',
         'Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Pseudoflavonifractor',
         'Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia',
         'Root;Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae',
         'Root;Bacteria;Bacteroidetes',
         'Root;Bacteria;Firmicutes;Clostridia',
         'Root;Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Alistipes')

for (otu in otus) {
  # make a short name (just the last rank) to use in the filename
  short_name = strsplit(otu, ';')[[1]] %>% tail(1)

  # make the filename where to save it
  fn = sprintf('timeseries-%s.pdf', short_name)

  # make and save the plot
  p = plot_timeseries_day14_ymax(otu) + ylab('relative abundance') + ggtitle(short_name)
  ggsave(fn, plot=p, useDingbats=F)
}
