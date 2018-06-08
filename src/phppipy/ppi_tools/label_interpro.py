#!/usr/bin/env python3
"""
Module to add GO annotations and re-map them to a desired depth.
"""

import collections

import numpy as np

from pathlib import Path


def create_uniprot2interpro_dict(filepath, uniprot_ac_list):
    """
    Creates a dictionary mapping UniProt ACs to InterPro domains.

    Parameters
    ----------
    filepath : the file path to a (filtered) InterPro map file.
    uniprot_ac_list : list
        A list of unique UniProt ACs present in the dataframe.

    Returns
    -------
    dict
        A dictionary mapping UniProt ACs to InterPro domains.

    """
    uniprot2interpro = collections.defaultdict(set)

    mapping_file = Path(filepath)

    with mapping_file.open() as mapping:
        for line in mapping:
            split_line = line.split('\t')
            if split_line[0] in uniprot_ac_list:
                # uniprot2interpro[split_line[0]] = split_line[1]
                uniprot2interpro[split_line[0]].add(split_line[1])
    return uniprot2interpro


def annotate_interpro(interaction_dataframe,
                      uniprot2interpro_dict,
                      columns=None):
    if not columns:
        columns = ['xref_A', 'xref_B']

    for i in columns:
        interaction_dataframe['interpro_' + i] = interaction_dataframe.apply(
            lambda x: uniprot2interpro_dict.get(x[i].split(':')[1], np.nan),
            axis=1)


#
# # Retrieve interpro domains
# def annotate_interpro(interaction_dataframe, xref_columns=list(('xref_A', 'xref_B'))):
#     entry_A = xref_columns[0]
#     entry_B = xref_columns[1]
#     entry_A_interpro = entry_A + '_interpro'
#     entry_B_interpro = entry_B + '_interpro'
#
#     # https://stackoverflow.com/questions/26977076/pandas-unique-values-multiple-columns
#
#     unique_ids = pd.Series(pd.unique(interaction_dataframe[[entry_A, entry_B]].values.ravel()))
#     unique_ids_split = '"' + unique_ids.str.split(':').str[1] + '"'
#     unique_ids_str = unique_ids_split.str.cat(sep=' ')
#
#     def joinit(iterable, delimiter):
#         it = iter(iterable)
#         yield next(it)
#         for x in it:
#             yield delimiter
#             yield x
#
#     unique_ids_split_extra = joinit(list(unique_ids_split), 'OR')
#     unique_ids_split_extra_str = ' '.join(unique_ids_split_extra)
#
#     import urllib.request;
#     urllib.parse
#
#     url = 'http://www.uniprot.org/uniprot/'
#
#     params = {
#         'format': 'tab',
#         'query': unique_ids_split_extra_str,
#         'columns' : 'id,database(interpro)'
#     }
#
#     data = urllib.parse.urlencode(params)
#     data = data.encode('ascii')
#
#     request = urllib.request.Request(url, data)
#
#     contact = "pieter.moris@uantwerpen.be"
#     request.add_header('User-Agent', 'Python %s' % contact)
#
#     with urllib.request.urlopen(request) as response:
#         mapping = response.read()
#
#     Path(savepath).parent.mkdir(parents=True, exist_ok=True)
#     Path(savepath).write_bytes(mapping)
#
#
#
#
#
#     interaction_dataframe[entry_A].str.extract('^.*:(\w*)-?', expand=False)
#     # interaction_dataframe[entry_A_interpro] =
#
#
#     interaction_dataframe['xref_A_interpro'] = interaction_dataframe['xref_A'].str.extract('^.*:(\w*)-?',
#                                                                                      expand=False).apply(
#         lambda x: gaf_dict.get(x, np.NaN))
#     interaction_dataframe['GO_B'] = interaction_dataframe['xref_B'].str.extract('^.*:(\w*)-?',
#                                                                                      expand=False).apply(
#         lambda x: gaf_dict.get(x, np.NaN))
