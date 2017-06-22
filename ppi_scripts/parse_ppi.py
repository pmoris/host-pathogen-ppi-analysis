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

s = read_mitab(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/hpidb2_March14_ 2017_mitab_plus.txt') '''


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
    return df


def read_mitab_virhost(filepath):
    file = Path(filepath)
    return pd.read_csv(file, sep='\t', encoding='ISO-8859-1',
                       names=['xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B', 'aliases_A',
                              'aliases_B', 'detection_method', 'author', 'publication', 'taxid_A', 'taxid_B',
                              'interaction_type', 'source_database_ids', 'interaction_identifiers', 'confidence_score'])


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


df_virhost = read_mitab_virhost(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/VirHostNet_January_2017.txt')

df_hpidb2 = read_mitab_hpidb2(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/hpidb2_March14_2017_mitab_plus.txt')

df_phisto = read_mitab_phisto(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/phisto_Jan19_2017.csv',
    r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/mi.obo')

# print(df_virhost.columns)
# print(df_hpidb2.columns)
# print(df_phisto.columns)
#
# print(df_hpidb2[0:1])
#
# df_hpidb2.columns.values[0] = 'xref_A'
#df.rename(columns={ df.columns[1]: "whatever" })

# print(df_hpidb2.valuecounts())

a_unique = df_hpidb2.taxid_A.unique()
b_unique = df_hpidb2.taxid_B.unique()

# print(a_unique)
print('length taxid_a',len(a_unique))
value = '9606'
if value in a_unique:
    print(value)
print(np.any([value in taxid for taxid in a_unique]))
# print(b_unique)
print('length taxid_b',len(b_unique))
print(np.any([value in taxid for taxid in b_unique]))


a_unique = df_virhost.taxid_A.unique()
b_unique = df_virhost.taxid_B.unique()

# print(a_unique)
print('length taxid_a',len(a_unique))
value = '9606'
if value in a_unique:
    print(value)
print(np.any([value in taxid for taxid in a_unique]))
# print(b_unique)
print('length taxid_b',len(b_unique))
print(np.any([value in taxid for taxid in b_unique]))

from taxid-data import *