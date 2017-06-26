#!/usr/bin/env python3

"""
Module to parse protein-protein interaction datasets in mitab format.
"""

import sys
import timeit

import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)

from pathlib import Path

'''
# Dask implementation
import dask.dataframe as dd
import dask.multiprocessing
def read_mitab_dask(filepath):
    file = Path(filepath)
    df = dd.read_csv(file, sep='\t', blocksize=1000000)
    df = df.compute(get=dask.multiprocessing.get)
    return df

s = read_mitab(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/'
               r'hpidb2_March14_2017_mitab_plus.txt')
'''

def read_mitab_hpidb2_manual(filepath):
    file = Path(filepath)
    with file.open(encoding="ISO-8859-1") as f:
        dict_list = []
        for line in f:
            l = line.strip().split('\t')
            keys = ['protein_xref_1', "protein_xref_2", "alternative_identifiers_1", "alternative_identifiers_2",
                    "protein_alias_1", "protein_alias_2", "detection_method", "author_name", "pmid", "protein_taxid_1",
                    "protein_taxid_2", "interaction_type", "source_database_id", "database_identifier", "confidence",
                    "protein_xref_1_unique", "protein_xref_2_unique", "protein_taxid_1_cat", "protein_taxid_2_cat",
                    "protein_taxid_1_name", "protein_taxid_2_name", "protein_seq1", "protein_seq2", "source_database",
                    "protein_xref_1_display_id", "protein_xerf_2_display_id"]
            dict_list.append({k: v for k, v in zip(keys, l)})
    df = pd.DataFrame(dict_list)
    return df


def read_mitab_pandas_chunk(filepath):
    file = Path(filepath)
    df = pd.read_csv(file, sep='\t', encoding="ISO-8859-1", iterator=True, chunksize=1000)
    df_concat = pd.concat([chunk for chunk in df])
    return df_concat


def read_mitab_hpidb2(filepath):
    # TODO: notification for regular or plus file...
    file = Path(filepath)
    df = pd.read_csv(file, sep='\t', encoding='ISO-8859-1', header=0)
    new_partial_columns = ['xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B', 'aliases_A',
                    'aliases_B', 'detection_method', 'author', 'publication', 'taxid_A', 'taxid_B',
                    'interaction_type', 'source_database_ids', 'interaction_identifiers', 'confidence_score']
    # re-name first few columns to same format used for VirHost data set.
    df.columns = new_partial_columns + df.columns.tolist()[len(new_partial_columns):]
    # rename_dict = dict(zip(df.columns.tolist()[:len(new_partial_columns)], new_partial_columns))
    # df.rename(columns=rename_dict, inplace=True)
    df.taxid_A = df.taxid_A.str.split('(').str.get(0)
    df.taxid_B = df.taxid_B.str.split('(').str[0]
    return df


def read_mitab_virhost(filepath):
    file = Path(filepath)
    df = pd.read_csv(file, sep='\t', encoding='ISO-8859-1',
                       names=['xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B', 'aliases_A',
                              'aliases_B', 'detection_method', 'author', 'publication', 'taxid_A', 'taxid_B',
                              'interaction_type', 'source_database_ids', 'interaction_identifiers', 'confidence_score'])
    df.detection_method = df.detection_method.str.replace('"', '')
    df.interaction_type = df.interaction_type.str.replace('"', '')
    df.source_database_ids = df.source_database_ids.str.replace('"', '')
    return df


def read_mitab_phisto(mitab_filepath, mi_filepath):
    file = Path(mitab_filepath)
    # display_id = entry name on uniprot
    # entry A is always pathogen partner in PHISTO
    df = pd.read_csv(file, sep=',', encoding='ISO-8859-1', header=0,
                     names=['pathogen', 'pathogen_taxid', 'xref_A', 'display_id_A',
                            'xref_B', 'display_id_B', 'detection_method', 'publication'])
    # Retrieve molecular interaction obo ontology
    reverse_mi_dict = {name: id for id, name in phisto_load_mi_ontology(mi_filepath).items()}
    # Fetch psi-mi id for given detection method in PHISTO interaction file,
    # otherwise, keep the original detection method name. Required because of e.g. "Other Methods".
    df.detection_method = df.detection_method.map(lambda x: reverse_mi_dict.get(x, x))
    return df


def phisto_load_mi_ontology(filepath):
    # Function to parse PSI-MI molecular interaction .obo ontology file into a dictionary.
    obo_file = Path(filepath)
    with obo_file.open() as obo:
        mi_dict = {}
        is_mi_term = False
        for line in obo:
            if line.startswith('[Term]'):
                is_mi_term = True
            else:
                if is_mi_term == True:
                    if line.startswith('id:'):
                        # Careful, id itself contains a colon as well, so split on ':\s'
                        mi_id = line.split(': ')[1].strip()
                    elif line.startswith('name:') and mi_id:
                        mi_dict[mi_id] = line.split(':')[1].strip()
                    elif not line.strip():
                        is_mi_term = False
                        mi_id = None
                    else:
                        continue
                else:
                    continue
    return mi_dict


df_virhost = read_mitab_virhost(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim'
                                r'/ppi_data/VirHostNet_January_2017.txt')

df_hpidb2 = read_mitab_hpidb2(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim'
                              r'/ppi_data/hpidb2_March14_2017_mitab_plus.txt')

df_phisto = read_mitab_phisto(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim'
                              r'/ppi_data/phisto_Jan19_2017.csv',
                              r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim'
                              r'/ppi_data/mi.obo')

print('Checking for presence of humans in taxid A and B columns')
print('taxid:9606' in df_hpidb2.taxid_A.values) # true
# values returns array, faster
# unique is array, but needs to uniqueify firt
# print('taxid:9606' in df_hpidb2.taxid_A.unique())
# print(df_hppid2.taxid_A.isin(['taxid:9606']))
# print(df_hpidb2.taxid_A.str.contains('9606').any())
# print(df_hpidb2.taxid_A.str.contains('9606').sum())
# print(np.any(['9606' in taxid for taxid in df_hpidb2.taxid_A.unique()]))
# print(np.any(['9606' in taxid for taxid in df_hpidb2.taxid_A.values]))
# print(len(df_hpidb2.taxid_A))
print('taxid:9606' in df_hpidb2.taxid_B.unique()) # false
print('taxid:9606' in df_virhost.taxid_A.unique()) # true
print('taxid:9606' in df_virhost.taxid_B.unique()) # true

print('Checking number of taxids in A and B')
print(df_hpidb2.groupby('taxid_A').size())
print(df_hpidb2.groupby('taxid_B').size())
print(df_virhost.groupby('taxid_A').size())
print(df_virhost.groupby('taxid_B').size())

import retrieve_taxids

taxdump_dir = Path(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/'
                   r'taxid_data/taxdump')
names = taxdump_dir / 'names.dmp'
# name2taxid, taxid2name = retrieve_taxids.parse_taxid_names(str(names))
# print(len(name2taxid))

print('Checking for intra-species interactions')
print(np.where(df_hpidb2.taxid_A == df_hpidb2.taxid_B)) # none
print(np.where(df_virhost.taxid_A == df_virhost.taxid_B)) # many
print(np.where(df_virhost.taxid_A == df_virhost.taxid_B)[0].size) # 28664 / 752190
print(df_virhost.loc[(df_virhost.taxid_A == df_virhost.taxid_B), ['taxid_A', 'taxid_B']].shape)
print(df_virhost.size)

print('Checking for human-human interactions')
print(df_hpidb2.loc[(df_hpidb2.taxid_A.str.contains('9606')) & (df_hpidb2.taxid_B.str.contains('9606')),
                     ['taxid_A', 'taxid_B']].shape) # none
print(df_virhost.loc[(df_virhost.taxid_A.str.contains('9606')) & (df_virhost.taxid_B.str.contains('9606')),
                     ['taxid_A', 'taxid_B']].shape) # 26732 / 28664


# import time
# B = []
# A = time.time()
#
# df_virhost['boundaries'] = np.where(df_virhost.taxid_A == df_virhost.taxid_B, 'intra-species', 'inter-species')
# B.append(time.time()-A)
# df_virhost['boundaries2'] = ['intra-species' if row['taxid_A'] == row['taxid_B']
#                              else 'inter-species' for index, row in df_virhost.iterrows()]
# B.append(time.time()-A)
# df_virhost['boundaries3'] = df_virhost.apply(
#                             lambda x: 'intra-species' if x['taxid_A'] == x['taxid_B'] else 'inter-species', axis=1)
# B.append(time.time()-A)
# df_virhost['boundaries4'] = ['intra-species' if row[0] == row[1] else 'inter-species' for row in
#                              zip(df_virhost['taxid_A'], df_virhost['taxid_B'])]
# B.append(time.time()-A)
# df_virhost.loc[(df_virhost.taxid_A == df_virhost.taxid_B), 'boundaries5'] = 'intra-species'
# df_virhost.loc[~(df_virhost.taxid_A == df_virhost.taxid_B), 'boundaries5'] = 'inter-species'
# B.append(time.time()-A)
# print(B)
# df_virhost.boundaries.equals(df_virhost.boundaries2)
df_virhost['inter-intra'] = np.where(df_virhost.taxid_A == df_virhost.taxid_B, 'intra-species', 'inter-species')

#TODO: check for duplicates between dataframes, taking into account aliases and alternative identifiers

a = df_hpidb2.loc[df_hpidb2.xref_A.isin(df_virhost.xref_A), ['xref_A']].xref_A.unique()
b = df_hpidb2.loc[df_hpidb2.xref_A.isin(df_virhost.xref_B), 'xref_A'].unique()
print(np.unique(np.append(a,b)).shape)
print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A) | df_hpidb2.xref_A.isin(df_virhost.xref_B)),
                    :].xref_A.unique().shape)
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A) | df_hpidb2.xref_B.isin(df_virhost.xref_B)),
                    :].xref_B.unique().shape)
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) | ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
                    :].xref_B.unique().shape)
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
                    :].xref_B.unique().shape)
print(df_hpidb2.loc[(df_hpidb2.xref_B.isin(df_virhost.xref_A)) | (df_hpidb2.xref_B.isin(df_virhost.xref_B)),
                    :].xref_B.unique().shape)

# entries in hpidb2 in either A or B that don't have any corresponding entry in virhost
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)),
                    :].xref_B.unique().shape)

print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)),
                    :].xref_A.unique().shape)

print(np.unique(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)),
                    ['xref_A', 'xref_B']]).size)

# virhost human in B, non-human in A
# print(df_virhost.loc[(~df_virhost.taxid_A.str.contains)])

print(df_virhost.head())

def annotate_origin(interaction_dataframe):
    # add viral: or host: string to each identifier
    return interaction_dataframe

# print(df_hpidb2.groupby('taxid_A').size())

#
# print(df_virhost[df_virhost.taxid_A == '9606'])
#
# # print(np.where(df_virhost.taxid_A == '9606' && df_virhost.taxid_B == '9606', 1, 0))
#
# print(df_virhost[df_virhost['taxid_A'] == 'taxid:9606'])
# print(df_hpidb2[df_hpidb2['taxid_A'] == 'taxid:9606'])
#
# print('\n\n\n\n\n\n HEAD')
#
# print(df_hpidb2.head())
#
# print(df_hpidb2.taxid_A)
# print('taxid:9606' in df_hpidb2.taxid_B.unique())
