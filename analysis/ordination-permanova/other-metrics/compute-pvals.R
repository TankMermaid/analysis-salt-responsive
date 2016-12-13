#!/usr/bin/env Rscript

library(dplyr)
library(vegan)

# read in the sample metadata
sample_data = read_tsv('../../../data/rdp_g.melt') %>%
  select(sample, diet, group, day) %>%
  distinct

adonis_df = function(matrix_fn, n_perms=99999) {
  dat = read_tsv(matrix_fn)

  # figure out which OTUs (i.e., rows and columns) to use
  good_samples = data_frame(sample=names(dat)) %>%
    left_join(sample_data, by='sample') %>%
    filter(day == 14)

  # keep only those rows and columns
  good_idx = names(dat) %in% good_samples$sample
  dat = dat[good_idx, good_idx] %>% as.dist

  adonis(dat ~ diet, data=good_samples, strata=good_samples$group, permutations=n_perms) %$%
    aov.tab %>%
    as_data_frame %>%
    mutate(fn=matrix_fn) %>%
    head(1)
}

fns = c('bray-curtis/bray-curtis.dat',
        'unifrac/unifrac-w.dat',
        'unifrac/unifrac-uw.dat')

data_frame(fn=fns) %>%
  rowwise %>%
  do(adonis_df(.$fn)) %>%
  select(fn, pvalue=`Pr(>F)`, SumsOfSqs, MeanSqs, F.Model, R2) %>%
  write_tsv('pvalues.tsv')
