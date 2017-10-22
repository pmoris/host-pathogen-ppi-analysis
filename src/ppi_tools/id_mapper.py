#!/usr/bin/env python3

"""
Module to create protein id mappings and apply them to a dataframe.
"""

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
        to_map = interaction_dataframe.loc[interaction_dataframe[col].str.contains(description), col]
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

    print('Created mapping files between UniProt ACs and', ids_str, 'in the data sets.\nStored in', savepath+'.')


def map2uniprot(interaction_dataframe, filepath=r'../../data/interim/mappings/', columns=list(('xref_A', 'xref_B'))):
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
    # $ tail -n +2 hpidb2_March14_2017_mitab.txt | cut -f1 | sed -r 's/(^.*):.*/\1/g' | sort -u
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
            # Store mapping into dictionary where non-UniProt keys map to lists of UniProt ACs.
            mapping_dict = {}
            for j in mapping_array:
                if not j[0] in mapping_dict:
                    mapping_dict[j[0]] = [j[1]]
                else:
                    mapping_dict[j[0]].append(j[1])
                    # TODO: not needed, always 1 to 1 mapping?

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
    print('Converted all identifiers to UniProt ACs.')
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
