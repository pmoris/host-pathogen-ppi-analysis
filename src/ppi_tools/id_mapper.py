#!/usr/bin/env python3

"""
Module to create protein id mappings and apply them to a dataframe.
"""
import urllib.request;
urllib.parse
import numpy as np
import pandas as pd

from pathlib import Path


def create_mapping_files(interaction_dataframe, from_id, description, savepath, columns):
    """Create mapping files between uniprot AC's and other identifiers.

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
        E.g. Entrez, ensemblgenomes, etc.
    savepath : string
        Filepath where to write the mapping files. (The default, defined by map2uniprot(), is a mapping directory in
        the ppi_data directory in the parent directory relative to where the script is called.)
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
        # note: use startswith to avoid mixing ensembl and embl ids...
        to_map = interaction_dataframe.loc[interaction_dataframe[col].str.startswith(description), col]
        ids = ids.append(to_map)
    ids = ids.reset_index(drop=True)
    ids = ids.str.split(':').str[1].unique()
    ids_str = ' '.join(np.char.mod('%s', ids))

    # http://www.uniprot.org/help/api_idmapping#id_mapping_python_example
    # https://docs.python.org/3/howto/urllib2.html
    # https://www.biostars.org/p/66904/

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

    savepath = Path(savepath)
    savepath.parent.mkdir(parents=True, exist_ok=True)
    savepath.write_bytes(mapping)

    print('Created mapping files between UniProt ACs and', ids_str,
          'in the data sets.\nStored in', str(savepath.resolve()) + '.')


def map2uniprot(interaction_dataframe, filepath=r'../../data/interim/mappings/', columns=None):
    """ Remap identifiers to UniProt AC's.

    Replaces non-UniProt AC's for the specified columns of a pandas DataFrame.
    Identifiers that could not be mapped are left unchanged.
    Only reviewed identifiers are considered.
    Identifiers that map to multiple UniProt AC's are concatenated into a long string separated by 'MULT_MAP'.
    The reason being: passing a list or tuple would upset any str.contains() or str.startswith()
    lookups from that point onwards (boolean array would contain NaNs).


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
    # $ cat <(tail -n +2 ../raw/ppi_data/hpidb2_March14_2017_mitab.txt | cut -f1) <(tail -n +2 ../raw/ppi_data/hpidb2_March14_2017_mitab.txt | cut -f2) | sed -r 's/(^.*?):.*/\1/g' | sort -u
    # hpidb: ddbj/embl/genbank, ensembl, ensemblgenomes, entrez, gene/locuslink, intact, refseq, uniprotkb
    # hpidb: ensembl, entrez gene/locuslink, intact, uniprotkb
    # virhostnet: refseq, uniprotkb

    if not columns:
        columns = ['xref_A', 'xref_B']

    identifiers = {'ddbj/embl/genbank:': 'EMBL_ID', 'ensembl:': 'ENSEMBL_ID', 'ensemblgenomes:': 'ENSEMBLGENOME_ID',
                   'entrez gene/locuslink:': 'P_ENTREZGENEID', 'refseq:': 'P_REFSEQ_AC'}
    #NOTE: UniProt REST API does not support intact EBI:identifiers.

    for db_tag, id_abbreviation in identifiers.items():
        # Define filepath (without forward slashes)
        path = Path(filepath) / (db_tag.replace('/', '').strip(':') + r'2uniprot.txt')

        # Re-run this if data set changes...
        # TODO: add options to force re-run
        if not Path(path).is_file():
            create_mapping_files(interaction_dataframe, id_abbreviation, db_tag, str(path), columns)

        # Read in file as array
        with path.open() as mapping:
            # TODO: possible to keep unreviewed identifiers? Will result in many more 1:n mappings...
            mapping_array = np.array([[line.split('\t')[i] for i in [0, 2]]
                                      for line in mapping if '\treviewed' in line])
            # Store mapping into dictionary where non-UniProt keys map to lists of UniProtACs.
            mapping_dict = {}
            for j in mapping_array:
                from_id = j[0]
                uniprot_acc = 'uniprotkb:' + j[1]
                if not from_id in mapping_dict:
                    mapping_dict[from_id] = [uniprot_acc]
                else:
                    # even when only selecting reviewed proteins, there are still 1:n mappings
                    # e.g. http://www.uniprot.org/uniprot/?query=yourlist:M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W&sort=yourlist:M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W&columns=yourlist%28M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W%29,isomap%28M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W%29,id,entry+name,reviewed,protein+names,genes,organism,length
                    mapping_dict[from_id].append(uniprot_acc)

                    # print(len(entrez2uniprot_array)) # 526
                    # print(np.unique(entrez2uniprot_array[:, 0], return_index=True)) # 510
                    # unique_indices = np.unique(entrez2uniprot_array[:, 0], return_index=True)[1]
                    # mask = np.ones(len(entrez2uniprot_array[:, 0]), dtype=bool)
                    # mask[unique_indices] = False
                    # print(entrez2uniprot_array[mask])

            def lambda_helper(row_entry):
                # retrieve id mapping if it exists, otherwise retain original id
                new_id = mapping_dict.get(row_entry.split(':')[1], row_entry)
                # all new mappings are in list format
                if type(new_id) == list:
                    # if multiple mappings, return as list
                    if len(new_id) > 1:
                        return 'MULT_MAP'.join(new_id)
                    else: # retrieve id as string from list
                        return new_id[0]
                # if no list is returned by dict.get(), original id is passed
                return new_id

            for col in columns:     # note: use startswith to avoid mixing ensembl and embl ids...
                mapping_selection = interaction_dataframe[col].str.startswith(db_tag)
                interaction_dataframe.loc[mapping_selection, col] = \
                    interaction_dataframe.loc[mapping_selection, col].apply(lambda x: lambda_helper(x))

    print('Converted all found identifiers to UniProt ACs.')
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

"""
def explode_mult(interaction_dataframe, columns=None):
    # 1) select relevant parts of data frame
    # 2) do list splitting into multiple rows
    # 3) append to ~original_selection DF
    # or make it so every identifier is in a list, even if it's a single one
"""

def remove_mult(interaction_dataframe, columns=None):
    """ Remove PPIs where either entry has multiple mappings to UniProtKB ids.
    These identifiers should be removed because otherwise this would artificially inflate their numbers during
    frequent item set mining.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein identifiers for PPIs.
    columns : list
        The names of the columns containing the identifiers that need to be removed.
        (The defaults are xref_A and xref_B).

    Returns
    -------
    None
        The subset of the DataFrame with PPIs whose partners have 1 unique mapping to UniProt ACs.
    """
    if not columns:
        columns = ['xref_A', 'xref_B']

    selection = (~interaction_dataframe[columns[0]].str.contains('MULT_MAP') &
                 ~interaction_dataframe[columns[1]].str.contains('MULT_MAP'))

    print('Omitted {} PPIs due to the existance of multiple mappings.'.format(np.sum(~selection)))
    interaction_dataframe = interaction_dataframe.loc[selection].apply(lambda x: x)
    interaction_dataframe = interaction_dataframe.reset_index(drop=True)
    return interaction_dataframe

# def lookup(interaction_dataframe, columns=None):
#     """ Remap identifiers to UniProt AC's.
#
#     Tries to replace non-UniProt AC's for the specified columns of a pandas DataFrame by looking in the unique
#     xref columns.
#     Only replaces entry name if a UniProt AC is found.
#
#     Parameters
#     ----------
#     interaction_dataframe : DataFrame
#         DataFrame containing protein identifiers that need to be remapped to UniProt Accesion Numbers.
#     columns : list
#         The names of the columns containing the identifiers that need to be remapped.
#         (Defaults to None, which uses ['xref_A', 'xref_B', 'protein_xref_1_unique', 'protein_xref_1_unique'].
#
#     Returns
#     -------
#     None
#         Modifies the supplied DataFrame in-place.
#
#     """
#     if not columns:
#         columns = ['xref_A', 'xref_B', 'protein_xref_1_unique', 'protein_xref_2_unique']
#
#     for col_xref, col_lookup in zip(columns[:2], columns[2:]):
#         mapping_selection = ( ~ interaction_dataframe[col_xref].str.contains('uniprot', case=False) &
#                                 interaction_dataframe[col_lookup].str.contains('uniprot', case=False))
#         interaction_dataframe.loc[mapping_selection, col_xref] = \
#             interaction_dataframe[mapping_selection, col_lookup].apply(lambda x: x)
