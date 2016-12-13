This folder has scripts and data for performing PERMANOVA tests using distance
metrics besides JSD.

- `bray-curtis/` has an R script that uses the OTU table to create a distance matrix
- `unifrac/` has a Python script and some more data

The R script `compute-pvals.R` uses the Bray-Curtis and UniFrac distance
matrices produced by running those two sets of scripts to in turn compute
p-values.
