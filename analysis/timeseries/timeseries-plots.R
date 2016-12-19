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
  
  point_offset = function(diet, day, size=0.2) {
    if_else(day %in% c(-1, 14),
            if_else(diet == 'N', -size, size),
            0.0)
  }
  boxplot_width = if_else(all_data$day %in% c(-1, 14), 0.5, 1.0)

  all_data %>%
    mutate(day=as.factor(day),
           diet=factor(diet, c('N', 'H')),
           point_x=as.numeric(day) + point_offset(diet, day)) %>% 
    ggplot(aes(x=factor(day), y=ra,  fill=diet)) +
      #stat_boxplot(data=hsd_data, outlier.shape=NA, color='orange', fill=NA) +
      #stat_boxplot(data=nd_data, outlier.shape=NA, color='black', fill=NA) +
      #geom_point(data=hsd_data, color='orange', position=position_jitter(h=0, w=0.1)) +
      #geom_point(data=nd_data, color='black', position=position_jitter(h=0, w=0.1)) +
      #stat_boxplot(outlier.shape=NA, fill=NA) +
      stat_boxplot(outlier.shape=NA, width=boxplot_width) +
      geom_point(aes(x=point_x), position=position_jitter(h=0, w=0.1)) +
      xlab('day') +
      ylab('relative abundance') +
      scale_color_manual(values=c('gray', 'black')) +
      scale_fill_manual(values=c('white', 'gray')) +
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

p = plot_timeseries(otu[1])
p
stop("it's cool")

for (otu in otus) {
  # make a short name (just the last rank) to use in the filename
  short_name = strsplit(otu, ';')[[1]] %>% tail(1)

  # make the filename where to save it
  fn = sprintf('timeseries-%s.pdf', short_name)

  # make and save the plot
  p = plot_timeseries_day14_ymax(otu) + ylab('relative abundance') + ggtitle(short_name)
  ggsave(fn, plot=p, useDingbats=F)
}
