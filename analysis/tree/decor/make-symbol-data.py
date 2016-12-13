#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

import argparse, math

def lines_to_dict(lines, drop=1):
    '''Convert a two-column, tab-separated file to a dictionary'''
    # drop the header
    for i in range(drop):
        next(lines)

    dat = {}
    for line in lines:
        otu, val = line.rstrip().split("\t")
        dat[otu] = float(val)

    return dat

def decor_lines(pval_dict, fc_dict, max_pval=0.05, min_opacity=0.0, target_size=50.0, min_log=3):
    if sorted(pval_dict.keys()) != sorted(fc_dict.keys()):
        raise RuntimeError("keys in two data sets are not matched")

    # get the taxa and data values all in the same order
    otus = sorted(pval_dict.keys())
    pvals = [pval_dict[x] for x in otus]
    fcs = [fc_dict[x] for x in otus]

    max_fc = max([abs(x) for x in fcs])

    # produce the header
    yield "DATASET_SYMBOL"
    yield "SEPARATOR SPACE"
    yield "DATASET_LABEL foo"
    yield "COLOR #ff0000"
    yield "MAXIMUM_SIZE {}".format(target_size)
    yield "DATA"

    for otu, pval, fc in zip(otus, pvals, fcs):
        # determine fill and opacity from pvalue
        if pval < max_pval:
            fill = 1
            opacity = 1.0 - (math.log10(pval) + min_log) / min_log
        else:
            fill = 0
            opacity = 1.0

        # determine coloring from direction of the fold change
        if fc < 0:
            rgb = '255,0,0'
        else:
            rgb = '0,0,255'

        color = 'rgba({},{})'.format(rgb, opacity)

        # determine size from the magnitude of the fold change
        size = target_size * (abs(fc) / max_fc)

        # produce an output line 
        yield " ".join([str(x) for x in ["'" + otu + "'", 2, size, color, fill, 1.0]])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='make circles with opacity and sizes to show enrichment')
    p.add_argument('pval', type=argparse.FileType('r'), help='p-value data file')
    p.add_argument('fc', type=argparse.FileType('r'), help='fold-change data file')
    args = p.parse_args()

    pval_dict = lines_to_dict(args.pval)
    fc_dict = lines_to_dict(args.fc)

    for line in decor_lines(pval_dict, fc_dict):
        print(line)
