'''
Script to perform a rough first filter on ppi datasets based on the
chosen pathogen ID.

- Generates all (pathogen) child taxon IDs by providing the taxdump directory
and desired parent taxid as arguments.

- Retrieves all associated host taxon IDs by searching through the PPI files.

- User should provide a directory containing .mitab files and optionally
a phi*.csv and mi.obo file.

All files matching the extensions in the given directory are imported,
in a non-recursive manner (i.e. sub-folders are not checked).

Output is saved to a directory named after the pathogen taxonid.
'''

import argparse
import numpy as np
import pandas as pd

from pathlib import Path

from phppipy.ppi_tools import ppi_import
from phppipy.dataprep import taxonid
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
    '-n',
    '--ncbi-taxdump',
    dest='taxdump',
    type=str,
    required=True,
    help='Taxdump directory with NCBI taxonomy files.')
parser.add_argument(
    '-t',
    '--taxonid',
    dest='taxonid',
    type=str,
    required=True,
    help='Pathogen taxon ID whose children will be retrieved.')
parser.add_argument(
    '-o',
    '--output',
    dest='output',
    type=str,
    required=True,
    help=
    'Output directory for saving child/host taxon ID and filtered PPI files.')
args = parser.parse_args()

# Retrieve child taxon IDs
try:
    taxdump_dir = Path(args.taxdump)
except IndexError:
    print('Incorrect path provided to taxdump directory.')
else:
    print('Parsing taxdump files...')
    names = taxdump_dir / 'names.dmp'
    name2taxid, taxid2name = taxonid.parse_taxid_names(names)
    nodes = taxdump_dir / 'nodes.dmp'
    taxid2parent, taxid2rank = taxonid.parse_taxid_nodes(nodes)
    parent2child = taxonid.create_parent2child_dict(taxid2parent)

    pathogen_taxid = args.taxonid
    print('Retrieving all child taxa of', pathogen_taxid)
    pathogen_child_ids = taxonid.get_children(pathogen_taxid, parent2child)

    # set output directory, suffixed with the pathogen taxid
    out_dir = Path(args.output) / pathogen_taxid
    out_path = out_dir / 'taxonid' / (pathogen_taxid + '-child-taxids.txt')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    taxonid.write_taxids(pathogen_child_ids, taxid2name, out_path)
    print('Saved output to', str(out_path), '\n')

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

# Filter on pathogen of choice
pathogen_child_ids_lookup = set(['taxid:' + i for i in pathogen_child_ids])
ppi_df = ppi_df.loc[(ppi_df.taxid_A.isin(pathogen_child_ids_lookup))
                    | (ppi_df.taxid_B.isin(pathogen_child_ids_lookup))]
ppi_df = ppi_df.reset_index(drop=True)

# Extract host taxids
# Remove ids already in pathogen id list, because the ppi datasets may contain
# intra-pathogen interactions too.
host_ids = {
    i
    for i in pd.unique(ppi_df[['taxid_A', 'taxid_B']].values.ravel('K'))
    if i not in pathogen_child_ids_lookup
}
# all_taxids = ppi_df['taxid_A'].append(ppi_df['taxid_B']).unique()
# host_ids = set(np.setdiff1d(all_taxids, list(pathogen_child_ids_lookup))) # don't forget to turn set back into list here
out_path = out_dir / 'taxonid' / 'interacting-taxids.txt'
taxonid.write_taxids([i.split(':')[1] for i in host_ids], taxid2name, out_path)
print('\nSaved host taxon IDs to', str(out_path))

# Print information about dataset

# Size information about dataset
print('\nSize of dataframe', ppi_df.shape)
print(ppi_df.groupby('origin').size())

## List associated taxon ids
print('\nThe following taxon IDs were found associated with {}:'.format(
    pathogen_taxid))
for i in host_ids:
    taxid = i.split(':')[1]
    count = ppi_df.loc[(ppi_df['taxid_A'] == i)
                       | (ppi_df['taxid_B'] == i)].shape[0]
    print(taxid, taxid2name.get(taxid, 'not found'), count)

## Check in which column pathogens occur
## phisto/hpidb store pathogens only in column B, others are mixed
for i in ppi_df.origin.unique():
    print('\nDataset {} contains pathogens in protein A?'.format(i),
          any(j in pathogen_child_ids_lookup
              for j in ppi_df.loc[ppi_df['origin'] == i, 'taxid_A']))
    print('Dataset {} contains pathogens in protein B?'.format(i),
          any(ppi_df.loc[ppi_df['origin'] == i, 'taxid_B'].isin(
              pathogen_child_ids_lookup)))

## Check the intra-pathogen interaction counts
print('\n Count of intra-pathogen interactions within species\n', ppi_df.loc[
    ppi_df.taxid_A.isin(pathogen_child_ids_lookup)
    & ppi_df.taxid_B.isin(pathogen_child_ids_lookup)].groupby('origin').size())

## Check the intra-host interaction counts for good measure
print('\nCount of intra-host interactions\n',
      ppi_df.loc[~(ppi_df.taxid_A.isin(pathogen_child_ids_lookup))
                 & ~(ppi_df.taxid_B.isin(pathogen_child_ids_lookup))].groupby(
                     'origin').size())

## Counts of inter-intra species and host-pathogen per data source
ppi_filter.annotate_inter_intra_species(ppi_df)
ppi_filter.annotate_inter_intra_pathogen(ppi_df, pathogen_child_ids_lookup)
ppi_df['taxid_name_A'] = ppi_df['taxid_A'].apply(
    lambda x: taxid2name.get(x.split(':')[1], 'missing name'))
ppi_df['taxid_name_B'] = ppi_df['taxid_B'].apply(
    lambda x: taxid2name.get(x.lstrip('taxid:'), 'missing name'))
# print(
#     'The following interactions are between different pathogen or host species'
# )
### print entries
# print(ppi_df.loc[
#     ppi_df['inter-intra-species'] != ppi_df['inter-intra-pathogen'], [
#         'taxid_A_name', 'taxid_B_name', 'inter-intra-pathogen',
#         'inter-intra-species', 'publication'
#     ]].drop_duplicates())
### Counts per dataset
print('\nCount of intra-pathogen interactions across species per dataset')
print(
    'intra-species != intra-path',
    ppi_df.loc[ppi_df['inter-intra-species'] != ppi_df['inter-intra-pathogen']]
    .groupby('origin').size())
# print('a=path & b=path & a!=b', ppi_df.loc[
#     (ppi_df.taxid_A != ppi_df.taxid_B)
#     & (ppi_df.taxid_A.isin(pathogen_child_ids_lookup)) &
#     (ppi_df.taxid_B.isin(pathogen_child_ids_lookup))].groupby('origin').size())
# print(ppi_df.loc[ (ppi_df['inter-intra-pathogen'] == 'intra') & (ppi_df.taxid_A !=
#                  ppi_df.taxid_B)].groupby('origin').size())

# save merged PPI dataset
out_path = out_dir / 'ppi_data' / (pathogen_taxid + '-ppi-merged.tsv')
out_path.parent.mkdir(parents=True, exist_ok=True)
ppi_df.to_csv(out_path, sep='\t', index=False, header=True)
print('\nSaved merged PPI dataset to {}'.format(out_path))

# import ipdb
# ipdb.set_trace()
