#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(vegan)

# read in the OTU table
otu_table = read_tsv('../data/rdp_g.melt') %>%
  # and compute relative abundances
  group_by(sample) %>%
  mutate(sample_counts=sum(counts), ra=counts/sample_counts) %>%
  ungroup() %>%
  filter(sample_counts >= 1000)

# cast this OTU table as a square matrix
otu_matrix = otu_table %>%
  select(otu, sample, ra) %>%
  spread(sample, ra) %>%
  select(-otu)

# also create just a list of the sample metadata
sample_metadata = otu_table %>%
  select(sample, day, diet, group) %>%
  filter(!duplicated(sample))

# create a matrix of JSD's using this square matrix
jsd = function(p, q) {
  m = 0.5 * (p + q)
  d1 = p * log2(p / m)
  d2 = q * log2(q / m)
  d1[is.na(d1)] = 0
  d2[is.na(d2)] = 0

  0.5 * (sum(d1) + sum(d2))
}

# create a table of JSDs
jsd_table = otu_table %$%
  # get all pairwise combinations of samples
  crossing(sample1=sample, sample2=sample) %>%
  # and compute the JSD for each pair of samples
  rowwise() %>%
  mutate(jsd=jsd(otu_matrix[[sample1]], otu_matrix[[sample2]])) %>%
  # bring in the metadata for each of the two samples
  left_join(sample_metadata, by=c('sample1'='sample')) %>%
  left_join(sample_metadata, by=c('sample2'='sample'), suffix=c('1', '2'))

# create the boxplot
p = jsd_table %>%
  # only single-count each comparison (i.e., keep only one entry
  # per distinct pairs of samples)
  filter(sample1 < sample2) %>%
  # keep only JSDs for day 14 samples
  filter(day1==14 & day2==14) %>%
  # evaluate the "type" of the comparison (e.g., between two ND samples?)
  mutate(type=factor(if_else(diet1 != diet2, 'ND vs. ND->HSD',
                     if_else(diet1 == 'H', 'ND->HSD vs. ND->HSD',
                     'ND vs. ND')),
                     levels=c('ND vs. ND', 'ND->HSD vs. ND->HSD', 'ND vs. ND->HSD'))) %>%
  ggplot(aes(y=jsd, x=type)) +
    geom_boxplot() +
    geom_point(position=position_jitter(h=0, w=0.1), cex=0.5) +
    ylim(0, 0.45) +
    ylab("Jensen-Shannon divergence") +
    xlab("") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1))

ggsave('jsd-group.pdf', width=89, height=80, units="mm", useDingbats=F)

# create a JSD distance matrix
jsd_matrix = jsd_table %>%
  select(sample1, sample2, jsd) %>%
  spread(sample2, jsd)

# ensure that the columns and rows have the same order
if (!all(jsd_matrix$sample1 == tail(colnames(jsd_matrix), -1))) {
  stop('jsd_matrix has mismatched row and column names')
}

jsd_metadata = jsd_matrix %>%
  select(sample=sample1) %>%
  left_join(sample_metadata, by='sample')

# turn jsd_matrix into a distance matrix
jsd_matrix %<>%
  select(-sample1) %>%
  as.dist

if (!all(jsd_metadata$sample == colnames(jsd_matrix))) {
  stop('jsd_metadata samples do not match matrix')
}

# get the subset of the matrix corresponding only to the
# day 14 comparisons
is_day14 = jsd_metadata$day == 14
jsd_metadata14 = jsd_metadata[is_day14, ]
jsd_matrix14 = jsd_matrix %>% as.matrix %>% .[is_day14, is_day14] %>% as.dist

# perform the PERMANOVA test and write its result to a file
res = adonis(jsd_matrix14 ~ diet, data=jsd_metadata14, strata=jsd_metadata14$group, permutations=99999)
capture.output(res, file='results.txt')

# perform the MDS ordination on the whole JSD matrix
ord = metaMDSiter(jsd_matrix)$points %>% as_data_frame
ord$type = jsd_metadata %$%
  if_else(diet == 'N' | day < 0, 'ND',
  if_else(diet == 'H' & day > 0 & day <= 14, 'ND->HSD',
  if_else(diet == 'H' & day > 14, 'ND->HSD->ND', 'other')))

# make the MDS plot
p = ord %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=type)) +
  theme_minimal() +
  scale_color_manual(values=c('black', 'orange', 'sky blue')) +
  coord_fixed() +
  theme(legend.position=(c(0.85, 0.85)),
        legend.background=element_rect(color='black'))

ggsave('mds.pdf', plot=p, useDingbats=F)
