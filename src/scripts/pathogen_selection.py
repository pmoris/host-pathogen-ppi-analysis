'''
Generate all (pathogen) child taxon IDs by providing the taxdump directory and
desired parent taxid as arguments.

Retrieves all associated host taxon IDs by searching through the PPI files.

Provide a directory containing .mitab files and optionally
a phi*.csv and mi.obo file.

All files matching the extensions are imported, in a non-recursive manner.
'''

import argparse
import numpy as np
import pandas as pd

from pathlib import Path

from phppipy.ppi_tools import import
from phppipy.dataprep import taxonid

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
    help='Output directory for saving child/host taxon ID files.')
args = parser.parse_args()

# Retrieve child taxon IDs
try:
    taxdump_dir = Path(args.taxdump)
except IndexError:
    print('Incorrect path provided.')
else:
    print('Parsing taxdump files...')
    names = taxdump_dir / 'names.dmp'
    name2taxid, taxid2name = taxonid.parse_taxid_names(names)
    nodes = taxdump_dir / 'nodes.dmp'
    taxid2parent, taxid2rank = taxonid.parse_taxid_nodes(nodes)
    parent2child = taxonid.create_parent2child_dict(taxid2parent)

    taxid = args.taxonid
    print('Retrieving all child taxa of', taxid)
    pathogen_child_ids = taxonid.get_children(taxid, parent2child)

    out_path = Path(args.output) / (taxid + '-childs.txt')
    taxonid.write_taxids(pathogen_child_ids, taxid2name, out_path)
    print('Saved output to', str(out_path), '\n')

# Import PPI files
input_dir = Path(args.input)
mitab_files = input_dir.glob('*.mitab')
ppi_df_list = [import.read_mi_tab(i) for i in mitab_files if i.is_file()]
phisto_files = input_dir.glob('phi*.csv')
mi_file = input_dir / 'mi.obo'
ppi_df_list.extend(
    [import.read_mitab_phisto(i, mi_file) for i in phisto_files])

# Merge PPI datasets
for i in ppi_df_list:
    i['origin'] = i.name
ppi_df = pd.concat(ppi_df_list, axis=0, join='outer', ignore_index=True)

# Filter on pathogens
pathogen_child_ids_lookup = set(['taxid:' + i for i in pathogen_child_ids])
ppi_df = ppi_df.loc[(ppi_df.taxid_A.isin(pathogen_child_ids_lookup))
                    | (ppi_df.taxid_B.isin(pathogen_child_ids_lookup))]
ppi_df = ppi_df.reset_index(drop=True)

# Size information about dataset
print('\nSize of dataframe', ppi_df.shape)
print(ppi_df.groupby('origin').size())

# Extract host taxids
host_ids = {
    i
    for i in pd.unique(ppi_df[['taxid_A', 'taxid_B']].values.ravel('K'))
    if i not in pathogen_child_ids_lookup
}
# all_taxids = ppi_df['taxid_A'].append(ppi_df['taxid_B']).unique()
# host_ids = set(np.setdiff1d(all_taxids, list(pathogen_child_ids_lookup))) # don't forget to turn set back into list here
out_path = Path(args.output) / 'taxid-hosts.txt'
taxonid.write_taxids([i.split(':')[1] for i in host_ids], taxid2name, out_path)
print('\nSaved host taxon IDs to', str(out_path))

# Print information about dataset

## List associated taxon ids
print('\nThe following taxon IDs were found associated with {}'.format(taxid))
for i in host_ids:
    taxid = i.split(':')[1]
    count = ppi_df.loc[(ppi_df['taxid_A'] == i)
                       | (ppi_df['taxid_B'] == i)].shape
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
print('\n Count of intra-pathogen interactions', ppi_df.loc[
    ppi_df.taxid_A.isin(pathogen_child_ids_lookup)
    & ppi_df.taxid_B.isin(pathogen_child_ids_lookup)].groupby('origin').size())

## Check the intra-host interaction counts for good measure
print('\nCount of intra-host interactions',
      ppi_df.loc[~(ppi_df.taxid_A.isin(pathogen_child_ids_lookup))
                 & ~(ppi_df.taxid_B.isin(pathogen_child_ids_lookup))].groupby(
                     'origin').size())

# save merged PPI dataset

# ppi_df.to_csv()


def annotate_inter_intra(interaction_dataframe):
    """Adds column to DataFrame specifying whether interaction is inter- or intra-species.

    The added column is named "inter-intra".

    Parameters
    ----------
    interaction_dataframe : DataFrame
        The pandas DataFrame should correspond to the PSI-MITAB format.

    Returns
    -------
    None
        Modifies DataFrame in-place by adding the "inter-intra" column.

    """
    interaction_dataframe['inter-intra'] = np.where(
        interaction_dataframe.taxid_A == interaction_dataframe.taxid_B,
        'intra-species', 'inter-species')


annotate_inter_intra(ppi_df)

ppi_df['inter-intra-general'] = np.where(
    ppi_df.taxid_A.isin(pathogen_child_ids_lookup) &
    ppi_df.taxid_B.isin(pathogen_child_ids_lookup), 'intra-species',
    'inter-species')

print(ppi_df[['inter-intra-general', 'inter-intra']])
print(all(ppi_df['inter-intra-general'] == ppi_df['inter-intra']))

ppi_df['taxid_A_name'] = ppi_df['taxid_A'].apply(
    lambda x: taxid2name.get(x.split(':')[1], 'missing'))
ppi_df['taxid_B_name'] = ppi_df['taxid_B'].apply(
    lambda x: taxid2name.get(x.lstrip('taxid:'), 'missing'))

print(ppi_df.loc[ppi_df['inter-intra-general'] != ppi_df['inter-intra'], [
    'taxid_A_name', 'taxid_B_name', 'inter-intra-general', 'inter-intra',
    'publication'
]].drop_duplicates())

import ipdb

ipdb.set_trace()