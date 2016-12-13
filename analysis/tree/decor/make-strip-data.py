#!/usr/bin/env python3

import argparse, re

def parse_template(template_lines):
    for line in template_lines:
        yield re.sub("#SHOW_INTERNAL 0", "SHOW_INTERNAL 1", line.rstrip())

def strip_data_lines(taxa, template_lines):
    for line in template_lines:
        yield line

    for taxon in taxa:
        yield " ".join([taxon, '#000000', taxon])


if __name__ == '__main__':
    p = argparse.ArgumentParser('Make a labeled strip for each taxon')
    p.add_argument('taxa', type=argparse.FileType('r'), help='list of taxa')
    p.add_argument('--template', '-t', type=argparse.FileType('r'), default=open('dataset_color_strip_template.txt'), help='template file')
    args = p.parse_args()

    taxa = ["'" + line.rstrip() + "'" for line in args.taxa]
    template_lines = parse_template(args.template)

    for line in strip_data_lines(taxa, template_lines):
        print(line)
