#!/usr/bin/env python3

"""
Script to filter relevant data out of the complete InterPro dataset (protein2ipr.dat).

Expects a new-line separated text file with UniProtACs.
"""

import argparse

from pathlib import Path

parser = argparse.ArgumentParser(
    description='Script to filter InterPro files (protein2ipr).',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    '-i',
    '--input',
    dest='input',
    type=str,
    required=True,
    help='Path to protein2ipr file.')
parser.add_argument(
    '-p',
    '--proteins',
    dest='proteins',
    type=str,
    required=True,
    help='Path to UniProtAC list.')
parser.add_argument(
    '-o',
    '--output',
    dest='output',
    type=str,
    required=True,
    help='Path to output file.')
args = parser.parse_args()

# import identifiers
infile_identifiers = Path(args.proteins)
with Path(infile_identifiers).open('r') as infile:
    identifiers = set([i.strip() for i in infile])

# import interpro file
protein2ipr_file = Path(args.input)

# create output directory
protein2ipr_reduced_file = Path(args.output)
protein2ipr_reduced_file.parent.mkdir(parents=True, exist_ok=True)

uniprot2interpro_dict = {}

with protein2ipr_file.open('r') as protein2ipr:
    with protein2ipr_reduced_file.open('w') as protein2ipr_reduced:
        outer_counter = 0
        inner_counter = 0
        for line in protein2ipr:
            if outer_counter % 1000000 == 0:
                print('Processed', outer_counter, 'lines')
            if line.split('\t')[0] in identifiers:
                inner_counter += 1
                protein2ipr_reduced.write(line)
                if inner_counter % 1000 == 0:
                    print('Written', inner_counter, 'lines.')
            outer_counter += 1
