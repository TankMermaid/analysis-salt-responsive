#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer) # to get colorRampPalette
library(grid) # to get unit

# make a distinctive color palette
qual_col_pals <- brewer.pal.info
qual_col_pals$names <- rownames(qual_col_pals)
qual_col_pals <- filter(qual_col_pals, category == "qual")
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, qual_col_pals$names))

# read in and prepare the OTU table for plotting
otu_table = read_tsv('../data/rdp_g.melt') %>%
  # compute relative abundances
  group_by(sample) %>% mutate(ra=counts/sum(counts)) %>% ungroup() %>%
  # exclude low-abundance OTUs from the plot
  filter(ra > 0.01) %>%
  # add an "animal id"
  mutate(animalid=interaction(diet, group, animal)) %>%
  # sort by names
  arrange(desc(otu))

# plot and save
taxa_plot = function(otu_table, output_fn, legend.size=0.3, legend.text=0.25, day.size=0.75) {
  p = otu_table %>%
    ggplot(aes(x=as.factor(day), y=ra, fill=otu)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=1.0) +
    geom_hline(yintercept=0.0) +
    facet_grid(. ~ animalid) +
    ylab("relative abundance") +
    xlab("day") +
    theme(legend.key.size=unit(legend.size, "cm"),
          legend.text=element_text(size=rel(legend.text)),
          axis.text.x=element_text(size=rel(day.size))) +
    scale_fill_manual(values=col_vector)

  ggsave(output_fn, p, width=50, units="cm")
}

taxa_plot(otu_table %>% filter(diet=="H", group==2), "hsd.pdf")
taxa_plot(otu_table %>% filter(diet=="N"), "nsd.pdf")
