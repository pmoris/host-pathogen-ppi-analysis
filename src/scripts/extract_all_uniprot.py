import argparse
import numpy as np
import pandas as pd
import sys

from pathlib import Path

from phppipy.ppi_tools import id_mapper
from phppipy.ppi_tools import ppi_import
from phppipy.ppi_tools import ppi_filter

parser = argparse.ArgumentParser(
    description='Script to extract taxon ids from ppi files.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    '-i',
    '--input',
    dest='input',
    type=str,
    required=True,
    help='Directory with PPI interaction files.')
parser.add_argument(
    '-m',
    '--mapping',
    dest='mapping',
    type=str,
    required=True,
    help='Full mapping file from EBI GOA.')
parser.add_argument(
    '-o',
    '--output',
    dest='output',
    type=str,
    required=True,
    help='Output path.')
args = parser.parse_args()

# Import PPI files
input_dir = Path(args.input)
mitab_files = input_dir.glob('*.mitab')
ppi_df_list = [ppi_import.read_mi_tab(i) for i in mitab_files if i.is_file()]
phisto_files = input_dir.glob('phi*.csv')
mi_file = input_dir / 'mi.obo'
ppi_df_list.extend(
    [ppi_import.read_mitab_phisto(i, mi_file) for i in phisto_files])

# Merge PPI datasets
for i in ppi_df_list:
    i['origin'] = i.name
ppi_df = pd.concat(ppi_df_list, axis=0, join='outer', ignore_index=True)

# remap to UniProt AC
id_mapper.check_unique_identifier(ppi_df)
out_mappings_dir = Path(args.output) / 'mapping'
full_mapping_file = Path(args.mapping)
try:
    out_mappings_dir.mkdir(parents=True, exist_ok=False)
    id_mapper.map2uniprot(ppi_df, out_mappings_dir, reviewed_only=True, full_mapping_file=full_mapping_file)
except FileExistsError:
    print(
        'Warning: supplied output directory already contains a "mapping" directory, aborting operation.'
    )
    sys.exit(1)

# remove multiple mappings
ppi_df = id_mapper.remove_mult(ppi_df)
# TODO: save these somewhere or deal with them in a way

# create unique identifier by combining xrefs
ppi_filter.unique_identifier(ppi_df)

# remove duplicates
ppi_df = ppi_df.drop_duplicates(subset=['xref_partners_sorted'], keep='first')
ppi_df = ppi_df.reset_index(drop=True)

# remove non-UniProt identifiers
ppi_df = ppi_df.loc[(ppi_df.xref_A.str.contains('uniprot'))
                    & (ppi_df.xref_B.str.contains('uniprot'))]
ppi_df = ppi_df.reset_index(drop=True)

# save protein list for filtering gaf/interpro files
all_identifiers = pd.Series(
    pd.unique(ppi_df[['xref_A',
                      'xref_B']].values.ravel('K'))).str.split(':').str.get(1)
out_identifiers = Path(args.output) / 'all_identifiers.txt'
out_identifiers.parent.mkdir(parents=True, exist_ok=True)

with out_identifiers.open('w') as out:
    for i in all_identifiers:
        out.write("{}\n".format(i))

print('Saved list of all UniProtACs to {}'.format(out_identifiers))
