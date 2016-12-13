#!/usr/bin/env Rscript

library(vegan)
library(readr)

# read in the OTU table
otu = read_tsv('../../../../data/rdp_g.counts') %>%
  # drop the OTU IDs
  select(-OTU_ID) %>%
  # perform the Bray-Curtis analysis
  t %>%
  vegdist(method='bray') %>%
  # transform it back into something easy to save
  as.matrix %>%
  as_data_frame %>%
  write_tsv('bray-curtis.dat')
