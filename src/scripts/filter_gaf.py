#!/usr/bin/env python3

"""
Script to filter relevant data out of the complete gene association file (goa_uniprot_all.gaf).

Expects a new-line separated text file with UniProtACs.
"""

import argparse

from pathlib import Path

parser = argparse.ArgumentParser(
    description='Script to filter gene association files.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    '-i',
    '--input',
    dest='input',
    type=str,
    required=True,
    help='Path to goa_uniprot_all.gaf file.')
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
gaf_file = Path(args.input)

# create output directory
gaf_reduced_file = Path(args.output)
gaf_reduced_file.parent.mkdir(parents=True, exist_ok=True)

uniprot2interpro_dict = {}

with gaf_file.open('r') as gaf:
    with gaf_reduced_file.open('w') as gaf_reduced:
        outer_counter = 0
        inner_counter = 0
        for line in gaf:
            if line.startswith('!'):
                continue
            if outer_counter % 10000000 == 0:
                print('Processed', outer_counter, 'lines')
            if line.split('\t')[1] in identifiers:
                inner_counter += 1
                gaf_reduced.write(line)
                if inner_counter % 1000 == 0:
                    print('Written', inner_counter, 'lines.')
            outer_counter += 1
