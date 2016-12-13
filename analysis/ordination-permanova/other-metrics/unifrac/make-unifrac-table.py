#!/usr/bin/env python2

import string, pandas as pd
from cogent.parse.tree import DndParser
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent.maths.unifrac.fast_tree import UniFracTreeNode

tr = string.maketrans(';', '|')

def table2dict(lines):
    '''Convert an OTU table into a nested dictionary of counts'''
    header_line = next(lines)
    header_fields = header_line.rstrip().split("\t")
    samples = header_fields[1:]

    dat = {}
    for line in lines:
        fields = line.rstrip().split("\t")
        otu = string.translate(fields[0], tr)
        counts = [int(x) for x in fields[1:]]
        dat[otu] = {s: c for s, c in zip(samples, counts)}

    return dat

with open('tree.newick') as f:
    raw_tree = f.read()

tree = DndParser(raw_tree, UniFracTreeNode)

with open('../../../../data/rdp_g.counts') as f:
    envs = table2dict(f)

# write the weighted and unweighted tables
for weighted, fn in [[True, 'unifrac-w.dat'], [False, 'unifrac-uw.dat']]:
    res = fast_unifrac(tree, envs, weighted=weighted)
    matrix, samples = res['distance_matrix']
    df = pd.DataFrame(data=matrix, index=samples, columns=samples)
    df.to_csv(fn, sep='\t', index=False)
