#!/usr/bin/env python3

"""
Module to parse protein-protein interaction datasets in mitab format.
"""

import os
import sys
import timeit

import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)

from pathlib import Path

sys.path.insert(0, '../go-enrichment-tool')
from go_enrichment_tool import obo_tools
from go_enrichment_tool import gaf_parser


def read_mitab_hpidb2(filepath):
    """Read HPIDB2 data set into pandas DataFrame.

    Renames headers to standard format used by this module.
    Can deal with both the standard and -plus data set provided by the HPIDB2.

    Parameters
    ----------
    filepath : string
        file path to HPIDB2 data set.

    Returns
    -------
    DataFrame
        DataFrame containing the HPIDB2 interaction data.

    """
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
    # name dataframe
    df.name = 'HPIDB2'
    return df


def read_mitab_virhost(filepath):
    """Read VirHost data set into pandas DataFrame.

    Renames headers to standard format used by this module.

    Parameters
    ----------
    filepath : string
        file path to VirHost data set.

    Returns
    -------
    DataFrame
        DataFrame containing the VirHost interaction data.

    """
    file = Path(filepath)
    df = pd.read_csv(file, sep='\t', encoding='ISO-8859-1',
                     names=['xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B', 'aliases_A',
                            'aliases_B', 'detection_method', 'author', 'publication', 'taxid_A', 'taxid_B',
                            'interaction_type', 'source_database_ids', 'interaction_identifiers', 'confidence_score'])
    df.detection_method = df.detection_method.str.replace('"', '')
    df.interaction_type = df.interaction_type.str.replace('"', '')
    df.source_database_ids = df.source_database_ids.str.replace('"', '')
    # name dataframe
    df.name = 'VirHostNet2'
    return df


def read_mitab_phisto(mitab_filepath, mi_filepath):
    """Read PHISTO data set into pandas DataFrame.

    Renames headers to standard format used by this module and modifies a number of entries to match the standard
    format.

    Parameters
    ----------
    mitab_filepath : string
        file path to PHISTO data set.
    mi_filepath : string
        file path to psi-mi .obo file.

    Returns
    -------
    DataFrame
        DataFrame containing the PHISTO interaction data.

    """
    file = Path(mitab_filepath)
    # display_id = entry name on uniprot
    # entry B is always pathogen partner in PHISTO
    df = pd.read_csv(file, sep=',', encoding='ISO-8859-1', header=0,
                     names=['pathogen', 'taxid_B', 'xref_B', 'display_id_B',
                            'xref_A', 'display_id_A', 'detection_method', 'publication'])
    # Retrieve molecular interaction obo ontology
    mi_dict = phisto_load_mi_ontology(mi_filepath)
    reverse_mi_dict = {name: id for id, name in mi_dict.items()}
    # Fetch psi-mi id for given detection method in PHISTO interaction file,
    # otherwise, keep the original detection method name. Required because of e.g. "Other Methods".
    df.detection_method = df.detection_method.map(lambda x: 'psi-mi:' + str(reverse_mi_dict.get(x, x)) +
                                                            '(' + str(x) + ')')
    # Append column name in front of entries
    df.publication = df.publication.map(lambda x: 'pubmed:' + str(x))
    df.xref_A = df.xref_A.map(lambda x: 'uniprotkb:' + str(x))
    df.xref_B = df.xref_B.map(lambda x: 'uniprotkb:' + str(x))
    # Add Human taxid_A column
    df['taxid_A'] = 'taxid:9606'
    # name dataframe
    df.name = 'PHISTO'
    return df


def phisto_load_mi_ontology(filepath):
    """ Read in PSI-MI molecular interaction ontology file and store in dictionary.

    Ontology file should be in .obo file. Can be retrieved from http://ontologies.berkeleybop.org/mi.obo


    Parameters
    ----------
    filepath : string
        The filepath to a PSI-MI ontology .obo file.

    Returns
    -------
    dictionary
        A dictionary mapping PSI-MI names to id's.

    """
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


def create_mapping_files(interaction_dataframe, from_id, description, filepath, columns):
    """ Create mapping files between uniprot AC's and other identifiers.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein identifiers that need to be remapped to UniProt Accesion Numbers.
    from_id : string
        The database abbreviations used by UniProt's mapping service, as described here:
        http://www.uniprot.org/help/api_idmapping
    description : string
        The naming convention used in the interaction data sets to describe the protein identifier source database.
        E.g. entrez, ensemblgenomes, etc.
    filepath : string
        Filepath where to write the mapping files. (The default, defined by map2uniprot(), is a ppi_data directory in
        the parent directory of where the script is called from, defined by .
    columns : list
        The names of the columns containing the identifiers that need to be remapped.
        (The defaults, defined by map2uniprot, are xref_A and xref_B).

    Returns
    -------
    None
        Writes mapping files to data directory.

    """
    # create space-separated Entrez gene string
    ids = pd.DataFrame()
    for col in columns:
        to_map = interaction_dataframe[col][interaction_dataframe[col].str.contains(description)]
        ids.append(to_map)
    ids = ids.reset_index(drop=True)
    ids = ids.str.split(':').str[1].unique()
    ids_str = ' '.join(np.char.mod('%s', ids))

    # http://www.uniprot.org/help/api_idmapping#id_mapping_python_example
    # https://docs.python.org/3/howto/urllib2.html
    # https://www.biostars.org/p/66904/

    import urllib.request;
    urllib.parse

    url = 'http://www.uniprot.org/uploadlists/'

    params = {
        'from': from_id,
        'to': 'ACC',
        'format': 'tab',
        'query': ids_str
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('ascii')

    request = urllib.request.Request(url, data)

    contact = "pieter.moris@uantwerpen.be"
    request.add_header('User-Agent', 'Python %s' % contact)

    with urllib.request.urlopen(request) as response:
        mapping = response.read()

    Path(filepath + description + r'2uniprot.txt').write_bytes(mapping)


def map2uniprot(interaction_dataframe, filepath=r'../ppi_data/', columns=['xref_A', 'xref_B']):
    """ Remap identifiers to UniProt AC's.

    Replaces non-UniProt AC's for the specified columns of a pandas DataFrame.
    Identifiers that could not be mapped are left unchanged.
    Identifiers that map to multiple UniProt AC's are concatenated into a tuple.
    Only reviewed identifiers are considered.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein identifiers that need to be remapped to UniProt Accesion Numbers.
    filepath : string
        Filepath where to write the mapping files. (The default is a ppi_data directory in the parent directory
        of where the script is called from.
    columns : list
        The names of the columns containing the identifiers that need to be remapped.
        (The defaults are xref_A and xref_B).

    Returns
    -------
    None
        Modifies the supplied DataFrame in-place.

    """
    #TODO: automatically retrieve descriptions and identifiers from file.
    for i in [['P_ENTREZGENEID', 'entrez'], ['EMBL_ID', 'embl'], ['ENSEMBLGENOME_ID', 'ensemblgenomes'],
              ['P_REFSEQ_AC', 'refseq']]:
        # Re-run this if data set changes...
        if not os.path.isfile(r'../ppi_data/' + i[1] + r'2uniprot.txt'):
            create_mapping_files(interaction_dataframe, i[0], i[1], filepath, columns)
        # Read in file as array
        with Path(r'../ppi_data/' + i[1] + r'2uniprot.txt').open() as mapping:
            # Skip unreviewed identifiers
            mapping_array = np.array([[line.split('\t')[i] for i in [0, 2, 3]]
                                      for line in mapping if '\treviewed' in line])
            # Store mapping into dictionary where non-uniprot keys map to lists of uniprot AC's.
            mapping_dict = {}
            for i in mapping_array:
                if not i[0] in mapping_dict:
                    mapping_dict[i[0]] = [i[1]]
                else:
                    mapping_dict[i[0]].append(i[1])

                    # print(len(entrez2uniprot_array)) # 526
                    # print(np.unique(entrez2uniprot_array[:, 0], return_index=True)) # 510
                    # unique_indices = np.unique(entrez2uniprot_array[:, 0], return_index=True)[1]
                    # mask = np.ones(len(entrez2uniprot_array[:, 0]), dtype=bool)
                    # mask[unique_indices] = False
                    # print(entrez2uniprot_array[mask])

            for col in columns:
                mapping_selection = (interaction_dataframe[col].str.contains(i[0])) | \
                                    (interaction_dataframe[col].str.contains('entrez'))

                interaction_dataframe.loc[mapping_selection, col] = interaction_dataframe.loc[
                    mapping_selection, col].apply(lambda x: tuple(mapping_dict[x.split(':')[1]]) if x in mapping_dict
                                                                                                 else x)


            # cat <(cut -f1 hpidb2_March14_2017_mitab.txt) <(cut -f2 hpidb2_March14_2017_mitab.txt) | grep entrez | sort -u |
            # sed -rn 's/^.*:(.*)/\1/p' | wc -l
            # 3044 entrez genes
            # entrez = df_hpidb2.loc[(df_hpidb2.xref_A.str.contains('entrez')) | (df_hpidb2.xref_B.str.contains('entrez'))]
            # entrez = entrez.xref_A.append(entrez.xref_B).unique()
            # print(entrez.shape)
            # cut -f1 test.txt | sort -u | wc -l
            # 523
            # grep -P '\treviewed' test.txt | cut -f1 | sort -u | wc -l
            # 510


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
    interaction_dataframe['inter-intra'] = np.where(interaction_dataframe.taxid_A == interaction_dataframe.taxid_B,
                                                    'intra-species', 'inter-species')
    # # Code comparing performance of various methods to add new column and fill in values.
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


def concat_interaction_datasets(list_of_datasets):
    """Concatenates a list of DataFrames containing molecular interaction information.

    Given a list of pandas DataFrames that conform to the PSI-MITAB format, concatenate them

    #and remove duplicated
    #interaction pairs. Duplicates are defined by a tuple of the "xref" entries (unique identifiers for the interactors).
    #Note that this will also omit duplicated pairs originating from different sources or found by other detection
    methods.

    Parameters
    ----------
    list_of_datasets : list
        A list of pandas DataFrames to be merged. The DataFrames should correspond to the PSI-MITAB format.

    Returns
    -------
    pd.DataFrame
        A new DataFrame built from the input DataFrames.

    """
    for i in list_of_datasets:
        i['origin'] = i.name

    df_concat = pd.concat(list_of_datasets, axis=0, ignore_index=True, join='outer')
    # df_concat = df_concat.drop_duplicates(subset=['xref_partners_sorted'])
    return df_concat



if __name__ == '__main__':
    # Read in PPI data sets
    df_virhost = read_mitab_virhost(r'../ppi_data/VirHostNet_January_2017.txt')

    df_hpidb2 = read_mitab_hpidb2(r'../ppi_data/hpidb2_March14_2017_mitab_plus.txt')

    df_phisto = read_mitab_phisto(r'../ppi_data/phisto_Jan19_2017.csv',
                                  r'../ppi_data/mi.obo')

    # TODO: no intact-EBI mapping?
    # Create files to map entrez gene id's to uniprot ac's
    # Map entrez gene id's to uniprot ac's

    map2uniprot(df_hpidb2)

# Concatenate the different sources
df_concat = concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto])

# Add unique identifier for interaction pairs
# df_concat['xref_partners_sorted'] = list(zip(df_concat.xref_A, df_concat.xref_B))
xref_partners_sorted_array = np.sort(np.stack((df_concat.xref_A, df_concat.xref_B), axis=1), axis=1)
df_concat['xref_partners_sorted'] = list(map(tuple, xref_partners_sorted_array))

# Label interactions as being within or between species.
annotate_inter_intra(df_concat)

# Remove duplicate interaction pairs (including different detection methods and publications)
df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted'])
df_concat_dedup = df_concat_dedup.reset_index(drop=True)

# Retrieve only Herpesviridae (taxid:10292), see retrieve_taxids.py script to generate child taxids
with Path(r'../taxid_data/child_taxids_of_10292.txt').open() as taxid_file:
    herpes_taxids = [str('taxid:' + line.split('|')[0]) for line in taxid_file]

# Filter herpes interactions
df_herpes = df_concat_dedup.loc[(df_concat_dedup.taxid_A.isin(herpes_taxids)) |
                                df_concat_dedup.taxid_B.isin(herpes_taxids)]
df_herpes = df_herpes.reset_index(drop=True)

print('\nNumber of inter versus intra interactions:\n')
print(df_herpes.groupby(['inter-intra', 'origin']).size())
print('\nNumber of taxid involved in intra-species interactions\n')
print(df_herpes.loc[df_herpes['inter-intra'] == 'intra-species'].groupby(['taxid_A', 'taxid_B']).size())  # no human

# Create gene ontology dataframe
obo_dict = obo_tools.importOBO('../go_data/go-basic.obo')
protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.split(':').str[1].unique())
gafDict = gaf_parser.importGAF('../go_data/gene_association_9606_10292.goa', protein_set)

print(df_herpes.loc[(df_herpes.xref_A.str.contains('entrez')) | (df_herpes.xref_B.str.contains('entrez'))].shape)

print(df_herpes.head())
df_herpes.to_csv('TEST.csv')

# Count missing values across columns
print('\nMissing values in each column:\n')
print(df_concat.isnull().sum(axis=0))

# Size of each data source
print('\nData source sizes:\n')
print(df_concat.groupby('origin').size())


# Check number of within and between interactions
def check_intra_species_interactions(interaction_dataframe):
    print('Checking for intra-species interactions')
    print(np.where(df_hpidb2.taxid_A == df_hpidb2.taxid_B))  # none
    print(np.where(df_virhost.taxid_A == df_virhost.taxid_B))  # many
    print(np.where(df_virhost.taxid_A == df_virhost.taxid_B)[0].size)  # 28664 / 752190
    print(df_virhost.loc[(df_virhost.taxid_A == df_virhost.taxid_B), ['taxid_A', 'taxid_B']].shape)


print('\nNumber of inter versus intra interactions:\n')
print(df_concat.groupby('inter-intra').size())
print(df_concat.groupby(['inter-intra', 'origin']).size())
print('\nNumber of taxid involved in intra-species interactions\n')
print(df_concat.loc[df_concat['inter-intra'] == 'intra-species'].groupby(['taxid_A', 'taxid_B']).size())
print('\nNumber of human-human interactions\n')
print(df_concat.loc[(df_concat['inter-intra'] == 'intra-species') & (df_concat['taxid_A'].str.contains('9606'))].shape)

# Count duplicates
print('\nNumber of duplicated interactions\n')
print(df_concat.duplicated(subset=['xref_partners_sorted']).shape)
print('\nNumber of unique interactions per data set\n')
print(df_concat.groupby('origin')['xref_partners_sorted'].nunique())


def drop_duplicates(interaction_dataframe, selection=[]):
    selection.extend(['xref_partners_sorted'])
    df = interaction_dataframe.drop_duplicates(subset=selection)
    return df


df_concat_dedup = drop_duplicates(df_concat,
                                  ['detection_method', 'publication'])  # 88645, 80855 for only detection method
print('\nSize of data set after removal of duplicates and before\n')
print(df_concat_dedup.shape)
print(df_concat.shape)

# TODO: check if duplicates from one data source differ in detection method, publication or something else
print('\nNumber of duplicate interaction pairs per data source \n')
print(df_concat.loc[df_concat.duplicated(subset=['xref_partners_sorted'], keep=False)].groupby('origin').size())


# TODO: when subsetting duplicates, also check taxid_A
# e.g. "Human herpesvirus 1 STRAIN KOS","10306","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
# "Human herpesvirus 1","10298","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
# TODO: differs in pathogen strain! Others are pubmed id based.


def check_host_location(interaction_dataframe):
    print('Checking for presence of humans in taxid A and B columns')
    print('taxid:9606' in df_hpidb2.taxid_A.values)  # true
    # values returns array, faster
    # unique is array, but needs to uniqueify firt
    # print('taxid:9606' in df_hpidb2.taxid_A.unique())
    # print(df_hppid2.taxid_A.isin(['taxid:9606']))
    # print(df_hpidb2.taxid_A.str.contains('9606').any())
    # print(df_hpidb2.taxid_A.str.contains('9606').sum())
    # print(np.any(['9606' in taxid for taxid in df_hpidb2.taxid_A.unique()]))
    # print(np.any(['9606' in taxid for taxid in df_hpidb2.taxid_A.values]))
    # print(len(df_hpidb2.taxid_A))
    print('taxid:9606' in df_hpidb2.taxid_B.unique())  # false
    print('taxid:9606' in df_virhost.taxid_A.unique())  # true
    print('taxid:9606' in df_virhost.taxid_B.unique())  # true


def count_taxids(interaction_dataframe):
    print('Checking number of taxids in A and B')
    print(df_hpidb2.groupby('taxid_A').size())  # contains non-human hosts
    print(df_hpidb2.groupby('taxid_B').size())
    print(df_virhost.groupby('taxid_A').size())
    print(df_virhost.groupby('taxid_B').size())


import retrieve_taxids

taxdump_dir = Path(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/'
                   r'taxid_data/taxdump')
names = taxdump_dir / 'names.dmp'


# name2taxid, taxid2name = retrieve_taxids.parse_taxid_names(str(names))
# print(len(name2taxid))



def check_intra_human_interactions(interaction_dataframe):
    print('Checking for human-human interactions')
    print(df_hpidb2.loc[(df_hpidb2.taxid_A.str.contains('9606')) & (df_hpidb2.taxid_B.str.contains('9606')),
                        ['taxid_A', 'taxid_B']].shape)  # none
    print(df_virhost.loc[(df_virhost.taxid_A.str.contains('9606')) & (df_virhost.taxid_B.str.contains('9606')),
                         ['taxid_A', 'taxid_B']].shape)  # 26732 / 28664


# TODO: check for duplicates between dataframes, taking into account aliases and alternative identifiers

print('HPIDB2 A proteins in VirHost A list')
a = df_hpidb2.loc[df_hpidb2.xref_A.isin(df_virhost.xref_A), ['xref_A']].xref_A.unique()
print('HPIDB2 A proteins in VirHost B list')
b = df_hpidb2.loc[df_hpidb2.xref_A.isin(df_virhost.xref_B), 'xref_A'].unique()
print('HPIDB2 A proteins in VirHost A or B list')
print(np.unique(np.append(a, b)).shape)
print('HPIDB2 A proteins in VirHost A or B list')
print(df_hpidb2.loc[(df_hpidb2.xref_A.isin(df_virhost.xref_A) | df_hpidb2.xref_A.isin(df_virhost.xref_B)),
      :].xref_A.unique().shape)
print('HPIDB2 A proteins unique count')
print(df_hpidb2.xref_A.unique().shape)
print('HPIDB2 A proteins not in VirHost A and not in B list (not in (A or B))')
print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A) | df_hpidb2.xref_A.isin(df_virhost.xref_B)),
      :].xref_A.unique().shape)
print('HPIDB2 B proteins not in VirHost A and not in B list')
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A) | df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins not in VirHost A or not in B list (')
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) | ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins not in VirHost A and not in B list (not in (A or B))')
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins in VirHost A or in B list')
print(df_hpidb2.loc[(df_hpidb2.xref_B.isin(df_virhost.xref_A)) | (df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins in VirHost A or in B list')
print(df_hpidb2.loc[(df_hpidb2.xref_B.isin(df_virhost.xref_A) | df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 A proteins not in VirHost A and not in B list OR B proteins not in VirHost A and not in B list')
# entries in hpidb2 in either A or B that don't have any corresponding entry in virhost
print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print(df_hpidb2.loc[(~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B))) |
                    (~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B))),
      :].xref_B.unique().shape)

print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_A.unique().shape)

print(
    np.unique(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)) |
                            ~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)),
                            ['xref_A', 'xref_B']]).size)


# virhost human in B, non-human in A
# print(df_virhost.loc[(~df_virhost.taxid_A.str.contains)])

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





#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# df1 = df_hpidb2.head().append(pd.DataFrame([['uniqueHpidb2_A', 'uniqueHpidb2_B'],
#                                             ['sharedAB_A', 'shared_AB_B'],
#                                             ['sharedA_A', 'sharedA_hpidb2B'],
#                                             ['sharedB_hpidb2A', 'sharedB_B'],
#                                             ['sharedReverse_A', 'sharedReverse_B'],
#                                             ['sharedA_virhostB', 'random'],
#                                             ['sharedB_B', 'sharedAB_A']],
#                                            columns=['xref_A', 'xref_B']), ignore_index=True)[df_hpidb2.columns.tolist()]
# df2 = df_virhost.head().append(pd.DataFrame([['uniqueVirhost_A', 'uniqueVirhost_B'],
#                                              ['sharedAB_A', 'shared_AB_B'],
#                                              ['sharedA_A', 'sharedA_virhostB'],
#                                              ['sharedB_virhostA', 'sharedB_B'],
#                                              ['sharedReverse_B', 'sharedReverse_A'],
#                                              ['sharedB_B', 'sharedB_virhostA']],
#                                             columns=['xref_A', 'xref_B']), ignore_index=True)[
#     df_virhost.columns.tolist()]
#
# # df2 = df2.append(pd.DataFrame([['uniqueVirhost_A', 'uniqueVirhost_B'],
# #                                ['sharedAB_A', 'shared_AB_B'],
# #                                ['sharedA_A', 'sharedA_virhostB'],
# #                                ['sharedB_virhostA', 'sharedB_B'],
# #                                ['sharedReverse_B', 'sharedReverse_A'],
# #                                ['sharedB_B', 'sharedB_virhostA']],
# #                                columns=['xref_A', 'xref_B']), ignore_index=True)[df_virhost.columns.tolist()]
# # df2.loc[12, 'extra'] = 'extra'
# print('testing\n\n\n\n\n')
# print(df1.loc[~(df1.xref_A.isin(df2.xref_A)) & ~(df1.xref_A.isin(df2.xref_B)) |
#               ~(df1.xref_B.isin(df2.xref_A)) & ~(df1.xref_B.isin(df2.xref_B)),
#       :])
#
# print('\n\ndf1\n\n\n')
# print(df1)
# print('\n\ndf2\n\n\n')
# print(df2)
# print('\n\nconcat\n\n\n')
#
# df1['xref_partners_sorted'] = list(zip(df1.xref_A, df1.xref_B))
# df2['xref_partners_sorted'] = list(zip(df2.xref_A, df2.xref_B))
#
# print(pd.concat([df1, df2], axis=0, ignore_index=True))
#
# df_merged = pd.merge(df1, df2, how='inner', on=['xref_partners_sorted'])
# print(df_merged)
#
# df_virhost['xref_partners_sorted'] = list(zip(df_virhost.xref_A, df_virhost.xref_B))
# df_hpidb2['xref_partners_sorted'] = list(zip(df_hpidb2.xref_A, df_hpidb2.xref_B))
# df_merged = pd.merge(df_virhost, df_hpidb2, how='inner', on=['xref_partners_sorted'])
# print(df_merged.xref_partners_sorted)
#
# print(len(df_merged['xref_partners_sorted'].unique()))
# print(len(df_merged['xref_partners_sorted']))
#
# print(len(df_hpidb2['xref_partners_sorted'].unique()))
# print(len(df_hpidb2['xref_partners_sorted']))
#
# print(len(df_virhost['xref_partners_sorted'].unique()))
# print(len(df_virhost['xref_partners_sorted']))
#
# df_merged = pd.concat([df_hpidb2, df_virhost], axis=0, ignore_index=True, join='outer')
# print(len(df_merged['xref_partners_sorted'].unique()))
# print(len(df_merged['xref_partners_sorted']))
# df_merged = df_merged.drop_duplicates(subset=['xref_partners_sorted'])
#
# print(len(df_merged['xref_partners_sorted'].unique()))
# print(len(df_merged['xref_partners_sorted']))
