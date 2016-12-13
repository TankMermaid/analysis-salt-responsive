# Analytics for salt & microbiome

## File structure

`data/` has:

- `rdp_g.counts`: an OTU table showing the counts of the RDP genus-level OTUs
- `rdp_g.melt`: a "melted" version of the OTU table including sample data

`analysis/` has subfolders corresponding to different analyses shown in the
paper:

- `ordination-permanova` has the ordination plot and PERMANOVA analyses (Figures S3A, S3B).
- `enrichment` has the analysis about which taxa were enriched in the two treatments. It has one file `enrichment.R` that performs the analysis.
- `tree` has the phylogenetic tree and the scripts used to decorate the tree with the enrichment information (Figure S3D).
- `classifier` has the classifier analysis (Figures 1A, 1B).
- `timeseries` has the analysis about the dynamics of specific taxa (Figure 1C). It has one file `timeseries-plots.R` that uses the OTU table to generate timeseries plots for select OTUs.

## Author

- Scott Olesen (swo [at] alum [dot] mit [dot] edu)
