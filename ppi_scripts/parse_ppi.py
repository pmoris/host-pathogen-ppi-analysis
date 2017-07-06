#!/usr/bin/env python3

"""
Module to parse protein-protein interaction datasets in mitab format.
"""

import os
import sys

import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)

from pathlib import Path

# sys.path.insert(0, '../go-enrichment-tool')
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


def create_mapping_files(interaction_dataframe, from_id, description, savepath, columns):
    """ Create mapping files between uniprot AC's and other identifiers.

    Queries the UniProt mapping service to retrieve mappings for all non-UniProt identifiers found in the dataset.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein identifiers that need to be remapped to UniProt Accesion Numbers.
    from_id : string
        The database abbreviations used by UniProt's mapping service, as described here:
        http://www.uniprot.org/help/api_idmapping
    description : string
        The naming convention used in the interaction datasets to describe the protein identifier source database.
        E.g. entrez, ensemblgenomes, etc.
    savepath : string
        Filepath where to write the mapping files. (The default, defined by map2uniprot(), is a mapping directory in
        the ppi_data directory in the parent directory relative to where the script is called..
    columns : list
        The names of the columns containing the identifiers that need to be remapped.
        (The defaults, defined by map2uniprot, are xref_A and xref_B).

    Returns
    -------
    None
        Writes mapping files to data directory.

    """
    # create space-separated Entrez gene string
    ids = pd.Series()
    for col in columns:
        to_map = interaction_dataframe[col][interaction_dataframe[col].str.contains(description)]
        ids = ids.append(to_map)
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

    Path(savepath).parent.mkdir(parents=True, exist_ok=True)
    Path(savepath).write_bytes(mapping)


def map2uniprot(interaction_dataframe, filepath=r'../ppi_data/mappings/', columns=['xref_A', 'xref_B']):
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
        # Define filepath
        path = filepath + i[1] + r'2uniprot.txt'
        # Re-run this if data set changes...
        if not Path(path).is_file():
            create_mapping_files(interaction_dataframe, i[0], i[1], path, columns)
        # Read in file as array
        with Path(path).open() as mapping:
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
                mapping_selection = (interaction_dataframe[col].str.contains(i[1])) | \
                                    (interaction_dataframe[col].str.contains(i[1]))

                interaction_dataframe.loc[mapping_selection, col] = interaction_dataframe.loc[
                    mapping_selection, col].apply(lambda x: tuple(mapping_dict[x.split(':')[1]]) if x in mapping_dict
                                                                                                 else x)

            '''
            # cat <(cut -f1 hpidb2_March14_2017_mitab.txt) <(cut -f2 hpidb2_March14_2017_mitab.txt) | grep entrez | sort -u |
            # sed -rn 's/^.*:(.*)/\1/p' | wc -l
            # 3044 entrez genes
            entrez = df_hpidb2.loc[(df_hpidb2.xref_A.str.contains('entrez')) | (df_hpidb2.xref_B.str.contains('entrez'))]
            entrez = entrez.xref_A.append(entrez.xref_B).unique()
            print(entrez.shape)
            # cut -f1 test.txt | sort -u | wc -l
            # 523
            # grep -P '\treviewed' test.txt | cut -f1 | sort -u | wc -l
            # 510
            '''


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
    # import timeit

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


def reorder_pathogen_host_entries(interaction_dataframe, host_list = ['taxid:9606']):
    ''' Moves all pathogen entries to B columns and host entries to A columns.

    Selects all interaction entries where host entries occur in B columns instead of A
    and swaps the A and B columns.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        The pandas DataFrame containing the protein-protein interactions that need to be sorted.
    host_list : list
        List of host taxid's to look for in the B column.
        (Default value is ['taxid:9606'], e.g. human.

    Returns
    -------
    None
        Modifies the interaction DataFrame inplace.

    '''

    host_position_mask = interaction_dataframe['taxid_B'].isin(host_list)
    # Note: HPIDB2 always has hosts as partner A and PHISTO has human/pathogen labeled too.
    column_names = ['xref', 'taxid', 'aliases', 'alt_identifiers']
    columns_to_swap = [name + label for name in column_names for label in ['_A','_B']]
    columns_after_swap = [name + label for name in column_names for label in ['_B', '_A']]
    interaction_dataframe.loc[host_position_mask, columns_to_swap] = interaction_dataframe.loc[host_position_mask,
                                                                                               columns_after_swap].values
    # Slower alternative
    # lambdafunc = lambda x: pd.Series([x['xref_B'],x['xref_A']])
    # df_herpes.loc[host_position_mask, ['xref_A', 'xref_B']] = df_herpes.loc[host_position_mask].apply(lambdafunc, axis=1)


def annotate_GO(interaction_dataframe, gaf_dict):
        """ Adds Gene Ontology terms to interaction dataset.

        Creates separate columns with Gene Ontology terms for both interaction partner in the interaction DataFrame.
        Note that the UniProt AC's need to be formatted to remove potential PTM labels, e.g. P04295-PRO_0000115940.

        Parameters
        ----------
        interaction_dataframe : DataFrame
            The pandas DataFrame containing the protein-protein interactions to be labelled.

        gaf_dict : dictionary
            A dictionary mapping UniProt AC's to Gene Ontology terms.

        Returns
        -------
        None
            Modifies the interaction dataframe in-place by adding Gene Ontology columns.
        """
        # interaction_dataframe['xref_A_GO'] = interaction_dataframe['xref_A'].apply(lambda x:
        #                                                                            gaf_dict.get(x.split(':')[1]))
        # interaction_dataframe['xref_B_GO'] = interaction_dataframe['xref_B'].apply(lambda x:
        #                                                                            gaf_dict.get(x.split(':')[1]))
        interaction_dataframe['xref_A_GO'] = interaction_dataframe['xref_A'].str.extract('^.*:(\w*)-?',
                                                                                         expand=False).apply(lambda x:
                                                                                        gaf_dict.get(x, np.NaN))
        interaction_dataframe['xref_B_GO'] = interaction_dataframe['xref_B'].str.extract('^.*:(\w*)-?',
                                                                                         expand=False).apply(lambda x:
                                                                                        gaf_dict.get(x, np.NaN))


def add_label_to_GO(go_set, taxid):
    """Adds host/virus labels to GO terms.

    Should be called via an apply function that calls it for both partners' GO set and taxids
    for each interaction in a dataset.
    Note: currently all hosts are labelled identically.


    Parameters
    ----------
    go_set : set
        The set of GO terms associated with a specific protein partner, stored in the xref_A/B_GO columns after
        running the annotate_GO() function on an interaction dataframe.
    taxid : string
        The taxid of the protein partners, stored in the taxid_A/B columns of the interaction dataframe.
        E.g. taxid:9606

    Returns
    -------
    string
        A string of labeled GO terms separated by commas.

    """
    if pd.isnull(go_set):
        return go_set
    else:
        if taxid in host_taxids:
            label = 'host-'
        else:
            label = 'virus-'

        labeled_list = []
        for i in go_set:
            labeled_list.append(label+i)
        return ','.join(labeled_list)


if __name__ == '__main__':
    # Read in PPI datasets
    df_virhost = read_mitab_virhost(r'../ppi_data/VirHostNet_January_2017.txt')

    df_hpidb2 = read_mitab_hpidb2(r'../ppi_data/hpidb2_March14_2017_mitab_plus.txt')

    df_phisto = read_mitab_phisto(r'../ppi_data/phisto_Jan19_2017.csv',
                                  r'../ppi_data/mi.obo')

    # Concatenate the different sources
    df_concat = concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto])

    # TODO: no intact-EBI mapping?
    # Map entrez gene id's to uniprot ac's
    map2uniprot(df_concat)

    # Size of each data source
    print('\nData source sizes:\n')
    print(df_concat.groupby('origin').size())

    # Add unique identifier for interaction pairs
    # df_concat['xref_partners_sorted'] = list(zip(df_concat.xref_A, df_concat.xref_B))
    xref_partners_sorted_array = np.sort(np.stack((df_concat.xref_A, df_concat.xref_B), axis=1), axis=1)
    df_concat['xref_partners_sorted'] = list(map(tuple, xref_partners_sorted_array))
    df_concat['xref_partners_sorted'] = pd.Series(map(tuple, xref_partners_sorted_array))
    # slower alternative, returns series of lists
    # f_concat.apply(lambda x: sorted(x[['xref_A', 'xref_B']]), axis=1)
    # other options
    # https://stackoverflow.com/questions/40187874/python-pandas-two-columns-with-same-values-alphabetically-sorted-and-stored
    # use df[['col1','col2']].values to get array and use np.sort
    # or add .tolist() to get list and use sorted(i) for i in list
    # or sort in-place df.values.sort(axis=1)
    # pandas.DataFrame.sort_values

    # Count duplicates
    # TODO: check if duplicates from one data source differ in detection method, publication or something else
    print('\nNumber of duplicated interactions\n')
    print(df_concat.duplicated(subset=['xref_partners_sorted']).shape)
    print('\nNumber of unique interactions per data set\n')
    print(df_concat.groupby('origin')['xref_partners_sorted'].nunique())

    # Check non-uniprot AC's
    print('\nNumber of interactions without UniProt AC\n')
    print(df_concat.loc[(~df_concat.xref_A.str.contains('uniprot')) |
                        (~df_concat.xref_B.str.contains('uniprot'))].groupby('origin').size())

    # Label interactions as being within or between species.
    annotate_inter_intra(df_concat)

    # Remove duplicate interaction pairs (including different detection methods and publications)
    df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted'])
    df_concat_dedup = df_concat_dedup.reset_index(drop=True)

    # Retrieve only Herpesviridae (taxid:10292), see retrieve_taxids.py script to generate child taxids
    #TODO: import from retrieve_Taxids and create on the spot
    #TODO: then combine this code into a function to subset a given dataframe for a given taxid and its children
    #TODO: and a host list or default to all hosts
    #TODO: note that this filtering will generally result in intra-host interactions being omitted, while retaining intra-viral ones
    with Path(r'../taxid_data/child_taxids_of_10292.txt').open() as taxid_file:
        herpes_taxids = [str('taxid:' + line.split('|')[0]) for line in taxid_file]

    # Filter herpes interactions
    # This omits all human-human interactions, but not intra-virus interactions
    df_herpes = df_concat_dedup.loc[(df_concat_dedup.taxid_A.isin(herpes_taxids)) |
                                    df_concat_dedup.taxid_B.isin(herpes_taxids)]
    df_herpes = df_herpes.reset_index(drop=True)

    # Check how many non-uniprot interactions are left
    print('\nNumber of non-UniProt AC interactions for Herpes interactions\n')
    print(df_herpes.loc[(~df_herpes.xref_A.str.contains('uniprot')) |
                        (~df_herpes.xref_B.str.contains('uniprot'))].groupby('origin').size())

    # Filter on uniprot AC's
    df_herpes = df_herpes.loc[(df_herpes.xref_A.str.contains('uniprot')) & (df_herpes.xref_B.str.contains('uniprot'))]
    df_herpes = df_herpes.reset_index(drop=True)

    # TODO: when subsetting duplicates, also check taxid_A
    # e.g. "Human herpesvirus 1 STRAIN KOS","10306","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
    # "Human herpesvirus 1","10298","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
    # TODO: differs in pathogen strain! Others differ in their pubmed id

    # Print some information about the interactions
    print('\nNumber of inter versus intra interactions:\n')
    print(df_herpes.groupby(['inter-intra', 'origin']).size())
    print('\nNumber of taxid involved in intra-species interactions\n')
    print(df_herpes.loc[df_herpes['inter-intra'] == 'intra-species'].groupby(['taxid_A', 'taxid_B']).size())  # no human
    print('\nNumber of human-human interactions\n')
    print(df_herpes.loc[(df_herpes['inter-intra'] == 'intra-species') &
                        (df_herpes['taxid_A'].str.contains('9606'))].shape)
    print('\nNumber of interactions without UniProt AC\n')
    print(df_herpes.loc[(~df_herpes.xref_A.str.contains('uniprot')) |
                        (~df_herpes.xref_B.str.contains('uniprot'))].shape)

    # Check which hosts are present
    print('\nNumber of interactions per host\n')
    all_taxids = df_herpes['taxid_A'].append(df_herpes['taxid_B']).unique()
    host_taxids = list(np.setdiff1d(all_taxids, herpes_taxids))
    import retrieve_taxids
    taxid_names_path = Path(r'../taxid_data/taxdump/names.dmp')
    name2taxid, taxid2name = retrieve_taxids.parse_taxid_names(str(taxid_names_path))
    for i in host_taxids:
        taxid = i.split(':')[1]
        count = df_herpes['xref_partners_sorted'].loc[(df_herpes['taxid_A'] == i) | (df_herpes['taxid_B'] == i)].shape
        print(taxid, taxid2name[taxid], count)

    # Move all host partners in xref_B to xref_A (only an issue for VirHostNet)
    # Note, also move ALL other associated labels...
    #TODO: currently only swaps taxid and xref, nothing else.
    reorder_pathogen_host_entries(df_herpes, host_taxids)
    '''
    # https://stackoverflow.com/questions/25792619/what-is-correct-syntax-to-swap-column-values-for-selected-rows-in-a-pandas-data

    # print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_A.values)
    # print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_B.values)
    #
    # print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_A').size())
    # print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_B').size())
    #
    # print(df_herpes.groupby('taxid_A').size())
    # print(df_herpes.groupby('taxid_B').size())
    '''

    print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_A.values)
    print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_B.values)

    print('\n\n\n\n\nvirhosttaxidssizesgrouped\n\n\n\n')
    print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_A').size())
    print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_B').size())

    # Count missing values across columns
    print('\nMissing values in each column:\n')
    print(df_herpes.isnull().sum(axis=0))

    # Create Gene Ontology dictionaries
    #TODO: create quickGO query using host list + herpesviridae
    obo_dict = obo_tools.importOBO('../go_data/go-basic.obo')
    # protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.split(':').str[1].unique())
    protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.extract('^.*:(\w*)-?',
                                                                                         expand=False).unique())
    gaf_dict = gaf_parser.importGAF('../go_data/gene_association_hosts_10292.goa', protein_set)

    # Number of proteins without GO annotation
    #TODO: PTM processing id's e.g. P04295-PRO_0000115940
    print('\nNumber of proteins without GO annotation\n')
    not_annotated = [i for i in protein_set if i not in gaf_dict]
    print(len(not_annotated))

    # Add GO annotations
    annotate_GO(df_herpes, gaf_dict)

    # Add virus/host labels to GO
    lambda_go_label = lambda x: pd.Series([add_label_to_GO(x['xref_A_GO'], x['taxid_A']),
                                           add_label_to_GO(x['xref_B_GO'], x['taxid_B'])])
    df_herpes[['xref_A_GO', 'xref_B_GO']] = df_herpes.apply(lambda_go_label, axis=1)

    print(df_herpes.head())

    # TODO: map GO to fixed level


    # Save to transaction database
    # TODO: create a separate transaction base per virus type + only inter?
    df_output = df_herpes.loc[:, ['xref_partners_sorted', 'xref_A_GO', 'xref_B_GO']]
    df_output.to_csv(r'ppi_go_transactions.csv', sep='\t', index=False)

    for virus in np.sort(df_herpes['taxid_B'].unique()):
        print(taxid2name[virus.split(':')[1]])
    print('\n\n\n\n\n\n\n\n\n\n\n')
    for virus in np.sort(df_herpes['taxid_A'].unique()):
        print(taxid2name[virus.split(':')[1]])

    print('settings',list(np.setdiff1d(all_taxids, host_taxids)))

    for i in np.sort(np.setdiff1d(all_taxids, host_taxids)):
        print(taxid2name[i.split(':')[1]])

'''
df_herpes.loc[df_herpes['origin'] == 'VirHostNet2', 'xref_A'].apply(lambda x: x.split(':')[1])
versus
df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].apply(lambda x: x['xref_A'].split(':')[1], axis=1)

'''