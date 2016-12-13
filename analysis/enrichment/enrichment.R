#!/usr/bin/env Rscript

library(dplyr)

otu_table = read_tsv('../../data/rdp_g.melt') %>%
  # add pseudocounts
  mutate(pcounts=counts + 1) %>%
  # compute relative abundances (with and w/o pseudocounts)
  group_by(sample) %>%
  mutate(ra=counts/sum(counts), pra=pcounts/sum(pcounts)) %>%
  ungroup() %>%
  filter(day %in% c(-2, -1, 14))

# compute the fold changes
fc = otu_table %>%
  mutate(type=if_else(day==14 & diet=='H', 'HSD', 'ND')) %>%
  group_by(otu, type) %>%
  summarize(median_pra=median(pra)) %>%
  spread(type, median_pra) %>%
  mutate(fold_change = log(HSD/ND)) %>%
  select(otu, fold_change)

# if counts are 0 before and after, p-value should be 1.0
# otherwise, use t.test
ttest_pvalue = function(x, y) if_else(max(x, y) == 0 && unique(y) == 0, 1.0, t.test(x, y)$p.value)

pv = otu_table %>%
  mutate(type=if_else(day==14 & diet=='H', 'HSD', 'ND')) %>%
  group_by(otu) %>%
  summarize(raw_pvalue=ttest_pvalue(ra[which(type=='HSD')], ra[which(type=='ND')])) %>%
  mutate(pvalue=p.adjust(raw_pvalue, method='BH')) %>%
  select(otu, pvalue)

write_tsv(fc, 'data-fc.txt')
write_tsv(pv, 'data-pval.txt')
