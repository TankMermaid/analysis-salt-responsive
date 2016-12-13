- `taxa.txt` is a list of the taxa in the OTU table
- `fixrank2newick.py` is a Python 3 script that converts the list of taxa into a tree (Newick format)

The subfolder `decor/` has scripts used to create the decor for the tree.

- The two files `dataset_*_template.txt` are template data files from ITOL
- `make_strip_data.py` uses the taxa names and the template to make the data for the strip of taxa labels
- `make_symbol_data.py` uses the data produced in the enrichment analysis to make the sizes and colors of the circles for the tree
