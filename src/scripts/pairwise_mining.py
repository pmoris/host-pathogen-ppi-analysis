import argparse
import ast

import numpy as np
import pandas as pd

from pathlib import Path

# from goscripts import obo_tools
# from goscripts import gaf_parser

from phppipy.mining import pairwise

parser = argparse.ArgumentParser(
    description='Script to add annotations to PPI datasets.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    '-i',
    '--input',
    dest='input',
    type=str,
    required=True,
    help='Path to cleaned PPI dataset.')
parser.add_argument(
    '-o',
    '--output',
    dest='output',
    type=str,
    required=True,
    help='Path to output file.')
args = parser.parse_args()

# import ppi file
ppi_file = Path(args.input)
ppi_df = pd.read_csv(
    ppi_file,
    sep='\t',
    header=0)
    # converters={
    #     'GO_xref_A': ast.literal_eval,
    #     'GO_xref_B': ast.literal_eval,
    #     'interpro_xref_A': ast.literal_eval,
    #     'interpro_xref_B': ast.literal_eval
    # }
print('PPIs were imported from {}\n'.format(ppi_file))

# convert "string-ified" sets back to real sets
# NOTE: ast.literal_eval cannot operate on empty sets (or set())
for i in ['GO_xref_A', 'GO_xref_B', 'interpro_xref_A', 'interpro_xref_B']:
    ppi_df.loc[~ppi_df[i].isnull(), i] = ppi_df.loc[~ppi_df[i].isnull(), i].map(ast.literal_eval)

merged_annotations_A = pairwise.merge_annotations(
    ppi_df, ['GO_xref_A', 'interpro_xref_A'])
labeled_merged_annotations_A = pairwise.add_hp_label(merged_annotations_A, 'h')
ppi_df['annotations_A'] = labeled_merged_annotations_A

merged_annotations_B = pairwise.merge_annotations(
    ppi_df, ['GO_xref_B', 'interpro_xref_B'])
labeled_merged_annotations_B = pairwise.add_hp_label(merged_annotations_B, 'p')
ppi_df['annotations_B'] = labeled_merged_annotations_B

# ppi_df = ppi_df[0:10]

results = pairwise.find_pairwise(ppi_df, ['annotations_A', 'annotations_B'])
pairwise.multiple_testing_correction(
    results, columns=['chi2_p-value', 'G_p-value', 'fisher_p-value'])

# import ipdb; ipdb.set_trace()

out_path = Path(args.output)
out_path.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(out_path, sep='\t', index=True, header=True)
print('\nSaved pairwise association results to {}'.format(out_path))
