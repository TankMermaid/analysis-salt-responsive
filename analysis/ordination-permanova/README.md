This folder has scripts and data for performing PERMANOVA tests.

The R script `ordination-permanova.R` compute the JSD distance matrix, produces
the JSD-by-treatment-group boxplot, produces the MDS ordination plot, and runs
the PERMANOVA test using the JSD matrix.

The subfolder `other-metrics` contains scripts and data used to run the
PERMANOVA tests using other distance metrics.

- `other-metrics/bray-curtis/` has an R script that uses the OTU table to create a distance matrix
- `other-metrics/unifrac/` has a Python script and some more data

The R script `other-metrics/compute-pvals.R` uses the Bray-Curtis and UniFrac
distance matrices produced by running those two sets of scripts to in turn
compute p-values.
