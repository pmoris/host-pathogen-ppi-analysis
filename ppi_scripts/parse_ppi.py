#!/usr/bin/env python3

"""
Module to parse protein-protein interaction datasets in mitab format and convert them into transaction datasets.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from go_enrichment_tool import gaf_parser
from go_enrichment_tool import obo_tools
import retrieve_taxids

# for interaction console
# import sys
# sys.path.insert(0, '../go-enrichment-tool')
# sys.path.append("/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi_scripts/go_enrichment_tool")
# sys.path.append("/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi_scripts")

pd.set_option('display.max_columns', None)
pd.options.display.max_colwidth = 200


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
    df.taxid_B = df.taxid_B.map(lambda x: 'taxid:' + str(x))
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
                if is_mi_term:
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
    # create space-separated id string
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


def map2uniprot(interaction_dataframe, filepath=r'../ppi_data/mappings/', columns=list(('xref_A', 'xref_B'))):
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
    # TODO: automatically retrieve descriptions and identifiers from file.
    # bash one-liner to retrieve data sources
    # cut -f1 hpidb2_March14_2017_mitab.txt | sed -r 's/(^.*):.*/\1/g' | sort -u
    for i in [['P_ENTREZGENEID', 'entrez'], ['EMBL_ID', 'embl'], ['ENSEMBLGENOME_ID', 'ensemblgenomes'],
              ['P_REFSEQ_AC', 'refseq']]:
        # Define filepath
        path = filepath + i[1] + r'2uniprot.txt'
        # Re-run this if data set changes...
        # TODO: add options to force re-run
        if not Path(path).is_file():
            create_mapping_files(interaction_dataframe, i[0], i[1], path, columns)
        # Read in file as array
        with Path(path).open() as mapping:
            # Skip unreviewed identifiers
            mapping_array = np.array([[line.split('\t')[i] for i in [0, 2, 3]]
                                      for line in mapping if '\treviewed' in line])
            # Store mapping into dictionary where non-uniprot keys map to lists of uniprot AC's.
            mapping_dict = {}
            for j in mapping_array:
                if not j[0] in mapping_dict:
                    mapping_dict[j[0]] = [j[1]]
                else:
                    mapping_dict[j[0]].append(j[1])

                    # print(len(entrez2uniprot_array)) # 526
                    # print(np.unique(entrez2uniprot_array[:, 0], return_index=True)) # 510
                    # unique_indices = np.unique(entrez2uniprot_array[:, 0], return_index=True)[1]
                    # mask = np.ones(len(entrez2uniprot_array[:, 0]), dtype=bool)
                    # mask[unique_indices] = False
                    # print(entrez2uniprot_array[mask])

            for col in columns:
                mapping_selection = (interaction_dataframe[col].str.contains(i[1])) | \
                                    (interaction_dataframe[col].str.contains(i[1]))

                interaction_dataframe.loc[mapping_selection,
                                          col] = interaction_dataframe.loc[mapping_selection,
                                                                           col].apply(lambda x:
                                                                                      tuple(mapping_dict[x.split(':')[
                                                                                          1]]) if x in mapping_dict
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

    Given a list of pandas DataFrames that conform to the PSI-MITAB format and concatenate them.

    NOTE: THE ORDER IN WHICH THE DATASETS ARE COMBINED DETERMINES THE ORDER IN WHICH DUPLICATES WILL BE RETAINED.

    # and remove duplicated interaction pairs.
    # Duplicates are defined by a tuple of the "xref" entries (unique identifiers for the interactors).
    # Note that this will also omit duplicated pairs originating from different sources or found by other detection
    # methods.

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


def reorder_pathogen_host_entries(interaction_dataframe, host_list=list(('taxid:9606'))):
    ''' Moves all pathogen entries to B columns and host entries to A columns.

    Selects all interaction entries where host entries occur in B columns instead of A and swaps the A and B columns.
    Or vice versa, where pathogen occurs in column A.

    https://stackoverflow.com/questions/25792619/what-is-correct-syntax-to-swap-column-values-for-selected-rows-in-a-pandas-data

    Parameters
    ----------
    interaction_dataframe : DataFrame
        The pandas DataFrame containing the protein-protein interactions that need to be sorted.
    host_list : list
        List of host taxid's to look for in the B column.
        (Default value is ['taxid:9606'], e.g. human.)

    Returns
    -------
    None
        Modifies the interaction DataFrame inplace.

    '''

    host_position_mask = interaction_dataframe['taxid_B'].isin(host_list)
    # Note: HPIDB2 always has hosts as partner A and PHISTO has human/pathogen labeled too.
    column_names = ['xref', 'taxid', 'aliases', 'alt_identifiers']
    columns_to_swap = [name + label for name in column_names for label in ['_A', '_B']]
    columns_after_swap = [name + label for name in column_names for label in ['_B', '_A']]
    interaction_dataframe.loc[host_position_mask,columns_to_swap] = interaction_dataframe.loc[host_position_mask,
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
    interaction_dataframe['xref_A_GO'] = interaction_dataframe['xref_A'].str.extract('^.*:(\w*)-?', expand=False).apply(
        lambda x: gaf_dict.get(x, np.NaN))
    interaction_dataframe['xref_B_GO'] = interaction_dataframe['xref_B'].str.extract('^.*:(\w*)-?', expand=False).apply(
        lambda x: gaf_dict.get(x, np.NaN))


def label_host_pathogen(interaction_dataframe, pathogen_set, columns=list(('xref_A_GO', 'xref_B_GO')),
                               taxid_columns=list(('taxid_A', 'taxid_B'))):
    ''' Adds virus/host label to selected columns in DataFrame.

    The columns and taxid lists should of equal length and ordered in the same manner. E.g. for each position
    in the columns list, the corresponding position in the taxid list should repeat the associated taxid column.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        Pandas DataFrame that should be modified.
    pathogen_set : set
        A set of pathogen taxids of the format "taxid:#####".
    columns : list
        A list of columns whose entries should be labelled.
        (Default = ['xref_A_GO', 'xref_B_GO'])
    taxid_columns : list
        A list of taxid columns corresponding to the order of the columns parameter.
        (Default = ['taxid_A', 'taxid_B']
    Returns
    -------
    None
        Modifies the pandas DataFrame inplace.

    '''
    if not len(columns) == len(taxid_columns):
        raise ValueError("Lists of columns to label should match taxid columns.")

    else:
        for col, taxid in zip(columns, taxid_columns):
            lambda_label = lambda x: label_term(x[col], x[taxid], pathogen_set)
            interaction_dataframe[col] = interaction_dataframe.apply(lambda_label, axis=1)

def label_term(label_set, taxid, pathogen_set):
    """Adds host/virus labels to GO terms and returns string representation.

    Should be called via an apply function that calls it for both partners' GO set and taxids
    for each interaction in a dataset.
    Note: currently all hosts are labelled identically.


    Parameters
    ----------
    label_set : set
        The set of GO terms associated with a specific protein partner, stored in the xref_A/B_GO columns after
        running the annotate_GO() function on an interaction DataFrame.
    taxid : string
        The taxid of the protein partners, stored in the taxid_A/B columns of the interaction DataFrame.
        E.g. taxid:9606

    Returns
    -------
    string
        A string of labeled GO terms separated by commas.

    """
    if pd.isnull(label_set):
        return label_set
    else:
        if taxid in pathogen_set:
            label = 'virus-'
        else:
            label = 'host-'

        labeled_list = []
        for i in label_set:
            labeled_list.append(label + i)
        return ','.join(labeled_list)


if __name__ == '__main__':
    # Read in PPI datasets
    df_virhost = read_mitab_virhost(r'../ppi_data/VirHostNet_January_2017.txt')

    df_hpidb2 = read_mitab_hpidb2(r'../ppi_data/hpidb2_March14_2017_mitab_plus.txt')

    df_phisto = read_mitab_phisto(r'../ppi_data/phisto_Jan19_2017.csv',
                                  r'../ppi_data/mi.obo')

    # Concatenate the different sources
    df_concat = concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto])

    # Map entrez gene id's to uniprot ac's
    map2uniprot(df_concat)
    # TODO: no intact-EBI mapping...



    ##############
    ## PRINTOUT ##
    ##############
    # Size of each data source
    print('\nData source sizes:\n')
    print(df_concat.groupby('origin').size())

    # Add unique identifier for interaction pairs
    xref_partners_sorted_array = np.sort(np.stack((df_concat.xref_A, df_concat.xref_B), axis=1), axis=1)
    df_concat['xref_partners_sorted'] = list(map(tuple, xref_partners_sorted_array))
    df_concat['xref_partners_sorted'] = pd.Series(map(tuple, xref_partners_sorted_array))

    # sorting the following does not work:
    # df_concat['xref_partners_sorted'] = list(zip(df_concat.xref_A, df_concat.xref_B))

    # slower alternative, returns series of lists
    # f_concat.apply(lambda x: sorted(x[['xref_A', 'xref_B']]), axis=1)
    # other options
    # https://stackoverflow.com/questions/40187874/python-pandas-two-columns-with-same-values-alphabetically-sorted-and-stored
    # use df[['col1','col2']].values to get array and use np.sort
    # or add .tolist() to get list and use sorted(i) for i in list
    # or sort in-place df.values.sort(axis=1)
    # pandas.DataFrame.sort_values



    ##############
    ## PRINTOUT ##
    ##############
    # Count duplicates
    # TODO: check if duplicates from one data source differ in detection method, publication or something else
    print('\nNumber of duplicated interactions\n')
    print(df_concat.duplicated(subset=['xref_partners_sorted']).shape)
    print('\nNumber of unique interactions per data set\n')
    print(df_concat.groupby('origin')['xref_partners_sorted'].nunique())

    # TODO: bash script to check overlap: comm <(cut -f3 -d, phisto_Jan19_2017.csv | sed 's/"//g' | sort -u ) <(cut -f2 hpidb2_March14_2017_mitab.txt | sed s/uniprotkb://g | sort -u)

    # Check non-uniprot AC's
    print('\nNumber of interactions without UniProt AC\n')
    print(df_concat.loc[(~df_concat.xref_A.str.contains('uniprot')) |
                        (~df_concat.xref_B.str.contains('uniprot'))].groupby('origin').size())

    # Label interactions as being within or between species.
    annotate_inter_intra(df_concat)
    # TODO: filter on inter-species interactions?

    # Remove duplicate interaction pairs (including different detection methods and publications)
    # https://stackoverflow.com/a/41650846
    # https://stackoverflow.com/questions/33042777/
    # Note that this will result in the first dataset (e.g. hpidb2) having priority over the others.
    df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted', 'taxid_B', 'taxid_A'], keep='first')
    df_concat_dedup = df_concat_dedup.reset_index(drop=True)
    # When subsetting duplicates, also check taxid_A and taxid_B
    # e.g. "Human herpesvirus 1 STRAIN KOS","10306","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
    # "Human herpesvirus 1","10298","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
    # These differ in pathogen strain, but it's the same interaction.
    # Will come into play when comparing groups and different groups have same interaction
    # Which is not the case here, since both are human herpesvirus 1 strains
    # Similar issue for pubmed ID, method, etc.

    # Retrieve Herpesviridae (taxid:10292) list, see retrieve_taxids.py script to generate child taxids
    # TODO: import from retrieve_Taxids and create on the spot
    # TODO: then combine this code into a function to subset a given dataframe for a given taxid and its children
    # TODO: and a host list or default to all hosts
    # TODO: note that this filtering will generally result in intra-host interactions being omitted,
    # TODO: while retaining intra-viral ones
    with Path(r'../taxid_data/child_taxids_of_10292.txt').open() as taxid_file:
        herpes_taxids = [str('taxid:' + line.split('|')[0]) for line in taxid_file]

    # Filter herpes interactions
    # This omits all human-human interactions, but not intra-virus interactions
    df_herpes = df_concat_dedup.loc[(df_concat_dedup.taxid_A.isin(herpes_taxids)) |
                                    df_concat_dedup.taxid_B.isin(herpes_taxids)]
    df_herpes = df_herpes.reset_index(drop=True)

    ##############
    ## PRINTOUT ##
    ##############
    # Check how many non-uniprot interactions are left
    print('\nNumber of non-UniProt AC interactions for Herpes interactions\n')
    print(df_herpes.loc[(~df_herpes.xref_A.str.contains('uniprot')) |
                        (~df_herpes.xref_B.str.contains('uniprot'))].groupby('origin').size())

    # Filter on UniProt AC's due to GO mapping being UniProt-based
    df_herpes = df_herpes.loc[(df_herpes.xref_A.str.contains('uniprot')) & (df_herpes.xref_B.str.contains('uniprot'))]
    df_herpes = df_herpes.reset_index(drop=True)

    ##############
    ## PRINTOUT ##
    ##############
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

    taxid_names_path = Path(r'../taxid_data/taxdump/names.dmp')
    name2taxid, taxid2name = retrieve_taxids.parse_taxid_names(str(taxid_names_path))

    for i in host_taxids:
        taxid = i.split(':')[1]
        count = df_herpes['xref_partners_sorted'].loc[(df_herpes['taxid_A'] == i) | (df_herpes['taxid_B'] == i)].shape
        print(taxid, taxid2name[taxid], count)

    # Create taxid dictionaries
    taxid_nodes_path = Path(r'../taxid_data/taxdump/nodes.dmp')
    taxid2parent, taxid2rank = retrieve_taxids.parse_taxid_nodes(str(taxid_nodes_path))
    parent2child = retrieve_taxids.create_parent2child_dict(taxid2parent)

    # Move all host partners in xref_B to xref_A (only an issue for VirHostNet)
    # Note, also move ALL other associated labels...
    reorder_pathogen_host_entries(df_herpes, host_taxids)

    ##############
    ## PRINTOUT ##
    ##############
    '''
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
    # TODO: create quickGO query using host list + herpesviridae
    obo_dict = obo_tools.importOBO('../go_data/go-basic.obo')
    # protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.split(':').str[1].unique())
    # Needs special extract pattern because of PTM processing id's e.g. P04295-PRO_0000115940
    protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.extract('^.*:(\w*)-?',
                                                                                               expand=False).unique())
    gaf_dict = gaf_parser.importGAF('../go_data/gene_association_hosts_10292.goa', protein_set)

    ##############
    ## PRINTOUT ##
    ##############
    # Number of proteins without GO annotation
    print('\nNumber of proteins without GO annotation\n')
    not_annotated = [i for i in protein_set if i not in gaf_dict]
    print(len(not_annotated))

    # Add GO annotations
    annotate_GO(df_herpes, gaf_dict)

    # Add virus/host labels to GO
    label_host_pathogen(df_herpes, herpes_taxids)



    # Map pathogen id to higher label
    pathogen_group_dict = {'bovine_ah1': '10320', 'bovine_hv1': '79889', 'epstein_barr': '10376', 'equid_av1': '10326',
                           'equid_gv2': '12657', 'gallid_av2': '10390', 'human_hsv1': '10298', 'saimiri_gv2': '10381',
                           'human_av2': '10310', 'human_av3': '10335', 'human_bv5': '10359', 'human_gv8': '37296',
                           'human_bv6A': '32603', 'human_bv6B': '32604', 'murid_bv1': '10366', 'murid_gv4': '33708',
                           'papiine_gv1': '106332', 'suid_av1': '10345'}

    for i, j in pathogen_group_dict.items():
        pathogen_group_dict[i] = [j] + retrieve_taxids.get_children(j, parent2child)


    # # ateline_gh3 = ['85618']
    # pathogen_group_dict['bovine_ah1'] = ['10320'] + retrieve_taxids.get_children('10320',parent2child)
    # # bovine_gh4 = ['10385'].extend(retrieve_taxids.get_children('10385',parent2child))
    # pathogen_group_dict['bovine_hv1'] = ['79889'] + retrieve_taxids.get_children('79889', parent2child)
    # # Eleph_hv1 = ['654902'] # actual parent taxa = 146015
    # pathogen_group_dict['epstein_barr'] = ['10376'] + retrieve_taxids.get_children('10376',parent2child)  # or human gammaherpesvirus 4
    # pathogen_group_dict['equid_av1'] = ['10326'] + retrieve_taxids.get_children('10326', parent2child)  # or equine hv1
    # pathogen_group_dict['equid_gv2'] = ['12657'] + retrieve_taxids.get_children('12657', parent2child)  # or equine hv2
    # pathogen_group_dict['gallid_av2'] = ['10390'] + retrieve_taxids.get_children('10390', parent2child)
    # pathogen_group_dict['human_hsv1'] = ['10298'] + retrieve_taxids.get_children('10298', parent2child)  # or human av1
    # pathogen_group_dict['saimiri_gv2'] = ['10381'] + retrieve_taxids.get_children('10381',parent2child)  # or herpesvirus saimiri
    # pathogen_group_dict['human_av2'] = ['10310'] + retrieve_taxids.get_children('10310', parent2child)
    # pathogen_group_dict['human_av3'] = ['10335'] + retrieve_taxids.get_children('10335',parent2child)  # or varicella-zoster
    # pathogen_group_dict['human_bv5'] = ['10359'] + retrieve_taxids.get_children('10359',parent2child)  # or human cytomegalovirus
    # pathogen_group_dict['human_gv8'] = ['37296'] + retrieve_taxids.get_children('37296',parent2child)  # also contains ape viruses?
    # pathogen_group_dict['human_bv6A'] = ['32603'] + retrieve_taxids.get_children('32603', parent2child)
    # pathogen_group_dict['human_bv6B'] = ['32604'] + retrieve_taxids.get_children('32604', parent2child)
    # pathogen_group_dict['murid_bv1'] = ['10366'] + retrieve_taxids.get_children('10366', parent2child)
    # pathogen_group_dict['murid_gv4'] = ['33708'] + retrieve_taxids.get_children('33708',parent2child)  # or murine hv 68
    # pathogen_group_dict['papiine_gv1'] = ['106332'] + retrieve_taxids.get_children('106332', parent2child)
    # pathogen_group_dict['suid_av1'] = ['10345'] + retrieve_taxids.get_children('10345',parent2child)  # or pseudorabies virus

    def pathogen_group_mapper(taxid, pathogen_group_dict):
        for pathogen_group, taxids in pathogen_group_dict.items():
            if taxid in taxids:
                return pathogen_group
        return np.NaN


    # lambda_pathogen_mapper = lambda x: [pathogen_group for pathogen_group, taxids in pathogen_group_dict.items()
    #                                     if taxid in taxids]

    df_herpes['pathogen_groups'] = df_herpes.apply(lambda x: pathogen_group_mapper(x['taxid_B'].split(':')[1],
                                                                                   pathogen_group_dict), axis=1)

    # TODO: map GO to fixed level


    # Save to transaction database
    # Note: Pathlib functionality is broken in Pandas 0.20
    output_directory = Path(r'../transaction_datasets/')
    output_directory.mkdir(exist_ok=True)

    # TODO: create a separate transaction base per virus type + only inter?
    df_output = df_herpes.loc[:, ['xref_partners_sorted', 'xref_A_GO', 'xref_B_GO']]
    df_output.reset_index(drop=True)
    df_output.to_csv(str(output_directory) + r'/ppi_go_transactions.csv', sep=',', index=False)

    # Separated datasets
    for i in df_herpes['pathogen_groups'].dropna().unique():
        df_output_grouped = df_herpes.loc[df_herpes['pathogen_groups'] == i,
                                          ['xref_partners_sorted', 'xref_A_GO', 'xref_B_GO']]
        df_output_grouped.reset_index(drop=True)
        df_output_grouped.to_csv(str(output_directory) + r'/' + str(i) + '.csv', sep=',', index=False)










    print(df_herpes.groupby('pathogen_groups').size())

    ##############
    ## PRINTOUT ##
    ##############
    for virus in np.sort(df_herpes['taxid_B'].unique()):
        print(taxid2name[virus.split(':')[1]])
    print('\n\n\n\n\n\n\n\n\n\n\n')
    for virus in np.sort(df_herpes['taxid_A'].unique()):
        print(taxid2name[virus.split(':')[1]])

        # TODO: column A contains some viruses due to intra-viral interactions. check if column A is being used to select
        # hosts at any point...
        # TODO: also: filter on inter?

    print('\n\n\n\n\n\nWAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAn\n\n\n\n\nn\n')

    print(sorted([taxid2name[i.split(':')[1]] for i in np.sort(np.setdiff1d(all_taxids, host_taxids))]))

    print(df_herpes.groupby('taxid_B').size())

    print('\n\n\n\n\n\nn\n\n\n\n\nn\n')
    print(np.setdiff1d(df_phisto['xref_A'].unique(), df_hpidb2['xref_A'].unique()).size)
    print('\n\n\n\n\n\nn\n\n\n\n\nn\n')
    print(np.setdiff1d(df_phisto['xref_B'].unique(), df_hpidb2['xref_B'].unique()).size)

    print('uniprotkb:Q91LX9' in df_phisto.values)
    print('uniprotkb:Q91LX9' in df_concat.values)
    print('uniprotkb:Q91LX9' in df_concat_dedup.values)
    print('uniprotkb:Q91LX9' in df_herpes.values)

    print(df_herpes.groupby('origin').size())
    df_herpes['pathogen_type'] = df_herpes['taxid_B'].apply(lambda x: taxid2name[x.split(':')[1]])
    print(df_herpes.groupby('pathogen_type').size())

    print(pd.Series(np.sort(df_herpes['pathogen_type'].unique())).apply(lambda x: [x, name2taxid[x]]))

    # TODO: subset intra interactions, uniquify taxid_A versus B and group by them. check for cross-species
    # TODO: ie A != B

    '''
pd.Series(np.sort(df_herpes['pathogen_type'].unique())).apply(lambda x: [x, name2taxid[x]])
0                                  [Ateline gammaherpesvirus 3, 85618]
1                                   [Bovine alphaherpesvirus 1, 10320]
2                                   [Bovine gammaherpesvirus 4, 10385]
3                                 [Bovine herpesvirus type 1.1, 79889]
4                 [Bovine herpesvirus type 1.1 (strain Cooper), 10323]
5                   [Bovine herpesvirus type 1.1 (strain Jura), 31518]
6                   [Bovine herpesvirus type 1.1 (strain P8-2), 10324]
7                    [Elephantid herpesvirus 1 (isolate Kiba), 654902]
8                             [Epstein-barr virus strain ag876, 82830]
9                                    [Equid alphaherpesvirus 1, 10326]
10                      [Equid herpesvirus type 2 strain 86/87, 82831]
11                    [Equine herpesvirus type 1 (strain AB4P), 31520]
12              [Equine herpesvirus type 1 (strain Kentucky A), 10329]
13                                  [Gallid alphaherpesvirus 2, 10390]
14                  [Herpes simplex virus (type 1 / strain 17), 10299]
15                   [Herpes simplex virus (type 1 / strain F), 10304]
16              [Herpes simplex virus (type 1 / strain Patton), 10308]
17                            [Herpesvirus saimiri (strain 11), 10383]
18                           [Herpesvirus saimiri (strain 488), 10384]
19                                   [Human alphaherpesvirus 1, 10298]
20                                   [Human alphaherpesvirus 2, 10310]
21                                   [Human alphaherpesvirus 3, 10335]
22                                    [Human betaherpesvirus 5, 10359]
23                                   [Human gammaherpesvirus 4, 10376]
24                                   [Human gammaherpesvirus 8, 37296]
25                             [Human herpesvirus 1 strain KOS, 10306]
26                             [Human herpesvirus 2 strain 333, 10313]
27                            [Human herpesvirus 2 strain HG52, 10315]
28                           [Human herpesvirus 3 strain Dumas, 10338]
29                    [Human herpesvirus 3 strain Oka vaccine, 341980]
30                           [Human herpesvirus 4 strain B95-8, 10377]
31                           [Human herpesvirus 5 strain AD169, 10360]
32                         [Human herpesvirus 5 strain Merlin, 295027]
33                           [Human herpesvirus 5 strain Towne, 10363]

34                            [Human herpesvirus 6 (strain GS), 10369]
35                   [Human herpesvirus 6 (strain Uganda-1102), 10370]
36                             [Human herpesvirus 6 strain Z29, 36351]

37                           [Human herpesvirus 8 strain GK18, 868565]
38                                [Human herpesvirus 8 type M, 435895]

39              [Marek's disease herpesvirus type 1 strain MD5, 10389]

40                                    [Murid betaherpesvirus 1, 10366]

41                                   [Murid gammaherpesvirus 4, 33708]

42                       [Murine cytomegalovirus (strain K181), 69156]
43                      [Murine cytomegalovirus (strain Smith), 10367]

44                                [Papiine gammaherpesvirus 1, 106332]

45                                    [Suid alphaherpesvirus 1, 10345]
46    [Suid herpesvirus 1 (strain Indiana-Funkhauser / Becker), 31523]
dtype: object
(8124, 35)
    
    
    
    pathogen_type
Ateline gammaherpesvirus 3                                    1
Bovine alphaherpesvirus 1                                    13
Bovine gammaherpesvirus 4                                     1
Bovine herpesvirus type 1.1                                   6
Bovine herpesvirus type 1.1 (strain Cooper)                   6
Bovine herpesvirus type 1.1 (strain Jura)                     5
Bovine herpesvirus type 1.1 (strain P8-2)                     3
Elephantid herpesvirus 1 (isolate Kiba)                       1
Epstein-barr virus strain ag876                            2160
Equid alphaherpesvirus 1                                      8
Equid herpesvirus type 2 strain 86/87                         4
Equine herpesvirus type 1 (strain AB4P)                       2
Equine herpesvirus type 1 (strain Kentucky A)                 6
Gallid alphaherpesvirus 2                                     4
Herpes simplex virus (type 1 / strain 17)                   925
Herpes simplex virus (type 1 / strain F)                      4
Herpes simplex virus (type 1 / strain Patton)                 1
Herpesvirus saimiri (strain 11)                               5
Herpesvirus saimiri (strain 488)                              4
Human alphaherpesvirus 1                                     34
Human alphaherpesvirus 2                                      3
Human alphaherpesvirus 3                                    133
Human betaherpesvirus 5                                       6
Human gammaherpesvirus 4                                    318
Human gammaherpesvirus 8                                    235
Human herpesvirus 1 strain KOS                                1
Human herpesvirus 2 strain 333                                8
Human herpesvirus 2 strain HG52                              14
Human herpesvirus 3 strain Dumas                              4
Human herpesvirus 3 strain Oka vaccine                       48
Human herpesvirus 4 strain B95-8                           2420
Human herpesvirus 5 strain AD169                             64
Human herpesvirus 5 strain Merlin                           117
Human herpesvirus 5 strain Towne                             21
Human herpesvirus 6 (strain GS)                               1
Human herpesvirus 6 (strain Uganda-1102)                      2
Human herpesvirus 6 strain Z29                                3
Human herpesvirus 8 strain GK18                             369
Human herpesvirus 8 type M                                  281
Marek's disease herpesvirus type 1 strain MD5                17
Murid betaherpesvirus 1                                     330
Murid gammaherpesvirus 4                                    432
Murine cytomegalovirus (strain K181)                         34
Murine cytomegalovirus (strain Smith)                         3
Papiine gammaherpesvirus 1                                    4
Suid alphaherpesvirus 1                                      59
Suid herpesvirus 1 (strain Indiana-Funkhauser / Becker)       4
dtype: int64
(8124, 35)
    
    '''

'''
df_herpes.loc[df_herpes['origin'] == 'VirHostNet2', 'xref_A'].apply(lambda x: x.split(':')[1])
versus
df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].apply(lambda x: x['xref_A'].split(':')[1], axis=1)

'''

# TODO: use mask objects for more general subsetting? http://growthintel.com/exploratory-data-analysis-in-pandas/
