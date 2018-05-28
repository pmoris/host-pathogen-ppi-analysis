#!/usr/bin/env python3
"""
Module to create protein id mappings and apply them to a dataframe.
"""

import urllib.request
import urllib.parse
import numpy as np
import pandas as pd

from pathlib import Path


def _create_mapping_files(array, from_id, description, savepath):
    """Create mapping files between uniprot ACs and other identifiers.

    Queries the UniProt mapping service to retrieve mappings for all
    non-UniProt identifiers found in the dataset.

    Parameters
    ----------
    array : array-like
        List of all identifiers that need to be remapped to
        UniProt Accession Numbers. Obtained from the pandas DataFrame and
        column list that is passed to the map2uniprot() function.
        E.g. [refseq:NP_001179, uniprotkb:Q16611, ...]
    from_id : string
        The database abbreviations used by UniProt's mapping service, as
        described here: http://www.uniprot.org/help/api_idmapping
    description : string
        The naming tag convention used in the interaction datasets to describe
        the protein identifier source database.
        E.g. entrez, ensemblgenomes, etc.
    savepath : string
        Filepath where to write the mapping file, passed along by the
        map2uniprot() function.

    Returns
    -------
    None
        Writes mapping files to savepath directory.

    """
    # create space-separated id string
    identifier_series = pd.Series(array)
    # note: use startswith to avoid mixing ensembl and embl ids...
    ids = identifier_series[identifier_series.str.startswith(description)]

    if ids.empty:
        print('No {} identifiers found in dataset.\n'.format(description))

    else:
        ids = ids.reset_index(drop=True)
        ids = ids.str.split(':').str[1].unique()
        ids_str = ' '.join(np.char.mod('%s', ids))

        # http://www.uniprot.org/help/api_idmapping#id_mapping_python_example
        # https://docs.python.org/3/howto/urllib2.html
        # https://www.biostars.org/p/66904/

        url = 'https://www.uniprot.org/uploadlists/'

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

        print('{} {} identifiers found in dataset.'.format(
            description, ids.size))

        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        savepath.write_bytes(mapping)

        print(
            'Created mapping file between UniProt ACs and {} in: {}.\n'.format(
                description.strip(':'), savepath.resolve()))


def _extract_identifiers(entries):
    """
    Extracts a unique list of identifiers from a group of identifiers that
    follow the PSI-MITAB format, e.g. ['intact:EBI-476586', dip:DIP-33878N'].

    Parameters
    ----------
    entries : array-like
        An array of identifiers.
    Returns
    -------
    list
        A list of the unique identifiers, except for UniProt AC.
    """

    identifier_list = list(
        pd.Series(entries).str.split(':').str.get(0).unique())
    identifier_list.remove('uniprotkb')
    return identifier_list


def _create_mapping_dict(path, reviewed_only):
    """
    Create a dictionary mapping a given type of identifiers to another one,
    based on a file created through the EBI mapping services (see
    _create_mapping_files()).
    
    Parameters
    ----------
    path : str
        The filepath to the mapping file.
    reviewed_only : bool
        Whether or not unreviewed mappings should be included.
    Returns
    -------
    dict
        A dictionary reflecting the identifier mapping in the file.
    """

    with Path(path).open() as mapping_file:

        # Store mapping into dictionary where non-UniProt keys map to lists of UniProtACs.
        reviewed_mapping_dict = {}
        unreviewed_mapping_dict = {}

        for line in mapping_file:
            # [0] contains the original identifier
            # [1] is an empty element in most mapping files
            # [2] contains the UniProt AC
            split_line = line.split('\t')
            from_id = split_line[0]
            uniprot_ac = 'uniprotkb:' + split_line[2]

            if '\treviewed' in line:
                if from_id not in reviewed_mapping_dict:
                    reviewed_mapping_dict[from_id] = [uniprot_ac]
                else:
                    # even when only selecting reviewed proteins, there are still 1:n mappings
                    # e.g. http://www.uniprot.org/uniprot/?query=yourlist:M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W&sort=yourlist:M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W&columns=yourlist%28M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W%29,isomap%28M2017111583C3DD8CE55183C76102DC5D3A26728BF70646W%29,id,entry+name,reviewed,protein+names,genes,organism,length
                    reviewed_mapping_dict[from_id].append(uniprot_ac)

            elif '\tunreviewed' in line:
                if from_id not in unreviewed_mapping_dict:
                    unreviewed_mapping_dict[from_id] = [uniprot_ac]
                else:
                    unreviewed_mapping_dict[from_id].append(uniprot_ac)

    # merge dictionaries if reviewed_only is false
    if not reviewed_only:
        # the order is extremely important, later keys will receive priority
        mapping_dict = {**unreviewed_mapping_dict, **reviewed_mapping_dict}
        # action = 'included in'
        # total = len(mapping_dict)
        #TODO check overlap!
    else:
        mapping_dict = reviewed_mapping_dict
        # action = 'excluded from'
        # total = len(reviewed_mapping_dict) + len(unreviewed_mapping_dict)
    '''
    TODO: the unreviewed/reviewed counts might be off because some identifiers have both a reviewed and unreviewed mapping.
    2703399		F8REB8	F8REB8_HHV1	unreviewed	US2 (Virion protein US2)	US2 HHV1gp087	Human herpesvirus 1 (HHV-1) (Human herpes simplex virus 1)	291
    2703399		P06485	US02_HHV11	reviewed	Protein US2	US2	Human herpesvirus 1 (strain 17) (HHV-1) (Human herpes simplex virus 1)	291
    consequently, the size of the dictionaries is not representative of the total number of remapped identifiers!
    print('{} {} out of {} succesfully mapped identifiers were unreviewed and {} the dataset.\n'.format(db_tag, len(unreviewed_mapping_dict), total, action))
    '''

    return mapping_dict


def map2uniprot(interaction_dataframe,
                filepath,
                columns=None,
                reviewed_only=True):
    """ Remap identifiers in specific columns of a DataFrame to UniProt ACs.

    Replaces non-UniProt ACs for the specified columns of a pandas DataFrame.
    Identifiers that could not be mapped are left unchanged.

    NOTE: Only reviewed identifiers are used by default.
    NOTE: Identifiers that map to multiple UniProt ACs are concatenated into a
          long string separated by 'MULT_MAP'.
          The reason being: passing a list or tuple would upset any
          str.contains() or str.startswith() lookups from that point onwards
          (i.e. the boolean array would contain NaNs).

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein identifiers that need to be remapped to
        UniProt Accesion Numbers.
    filepath : string
        Filepath where to write the mapping files.
    columns : list
        The names of the columns containing the identifiers that need to be
        remapped.
        (The defaults are xref_A and xref_B).
    reviewed_only : boolean
        Specifies whether or not unreviewed mappings should be retained.
        (Default: True)

    Returns
    -------
    None
        Modifies the supplied DataFrame in-place.
    """
    # TODO: automatically retrieve descriptions and identifiers from file.
    # bash one-liner to retrieve data sources
    # $ cat <(tail -n +2 data/raw/ppi_data/*.mitab | cut -f1) <(tail -n +2 data/raw/ppi_data/*.mitab | cut -f2) | sed -r 's/(^[^:]*?):.*/\1/g' | sort -u | grep -v "|"
    # NOTE: some database entries contain multiple xrefs (e.g. dip:DIP-10051N|refseq:NP_285937|uniprotkb), hence the grep -v (technically not required since regex now matches
    # everything except colon until a colon, i.e. the first "dip" would be found, not the later refseq or uniprots in the same element).
    # chebi
    # ddbj/embl/genbank
    # ensembl
    # ensemblgenomes
    # entrez gene/locuslink
    # intact
    # refseq
    # uniprotkb
    # NOTE: UniProt REST API does not support intact EBI:identifiers.

    if not columns:
        columns = ['xref_A', 'xref_B']

    merged_columns = pd.unique(
        interaction_dataframe[columns].values.ravel('K'))
    present_identifiers = _extract_identifiers(merged_columns)

    identifier_dict = {
        'ddbj/embl/genbank:': 'EMBL_ID',
        'ensembl:': 'ENSEMBL_ID',
        'ensemblgenomes:': 'ENSEMBLGENOME_ID',
        'entrez gene/locuslink:': 'P_ENTREZGENEID',
        'refseq:': 'P_REFSEQ_AC',
        'dip:': 'DIP_ID'
    }

    for i in present_identifiers:
        if i + ':' not in identifier_dict:
            print(
                'WARNING: interaction dataset contains "{}" entries, which could not be remapped to UniProt AC (check all possible mappings at https://www.uniprot.org/help/api_idmapping .\n'.
                format(i))

    # for each non-uniprotac identifier defined in the dictionary:
    for db_tag, id_abbreviation in identifier_dict.items():

        # Define filepath (without forward slashes)
        path = Path(filepath) / (
            db_tag.replace('/', '-').strip(':') + r'2uniprot.txt')

        # create a mapping file using the uniprot mapping service
        # NOTE: does not create a file if the remapping is empty
        # Skip if file already exists TODO: add options to force re-run
        if not path.is_file():
            _create_mapping_files(merged_columns, id_abbreviation, db_tag,
                                  path)
        else:
            print(
                '{} mapping file already exists, not regenerating...\n'.format(
                    path.resolve()))

        # if a mapping file was created
        if Path(path).is_file():

            # create dictionaries for current identifier remapping file
            mapping_dict = _create_mapping_dict(path, reviewed_only)

            # create helper function to use in apply() acting on columns to be remapped
            def lambda_helper(row_entry, mapping_dict):
                # retrieve id mapping if it exists, otherwise retain original id
                original_id = row_entry.split(':')[1]
                new_id = mapping_dict.get(original_id, None)

                # if no list is returned by dict.get(), original id is passed
                if new_id:
                    # if a new id is found in the mapping_dict, return it
                    # mapping_dict contains lists, so extract or join for multiple mappings
                    if len(new_id) > 1:
                        new_id = 'MULT_MAP'.join(new_id)
                    else:  # retrieve id as string from list
                        new_id = new_id[0]
                else:
                    # return original id if no remapping is found
                    new_id = row_entry  # ROW ENTRY NOT original_id!!!
                return new_id

            # go through all IDs in the supplied columns and remap them if needed
            success_series = []
            total_series = []
            for col in columns:
                # note: use startswith to avoid mixing ensembl and embl ids...
                mapping_selection = interaction_dataframe[col].str.startswith(
                    db_tag)

                # keep track of total (unique) identifier count
                total_series.append(interaction_dataframe.loc[mapping_selection,col])

                interaction_dataframe.loc[mapping_selection, col] = \
                    interaction_dataframe.loc[mapping_selection, col].apply(func=lambda_helper, args=(mapping_dict,))

                # count how many (unique) identifiers were remapped
                success_series.append(interaction_dataframe.loc[mapping_selection & ~interaction_dataframe[col].str.startswith(db_tag), col])

            # find unique identifiers in total and successful remapped series
            success_count = pd.concat(success_series).unique().size
            total_count = pd.concat(total_series).unique().size

            print(
                '{} {} out of identifiers {} were succesfully remapped to UniProt accession numbers.\n'.
                format(db_tag, success_count, total_count))

    # report how many identifiers had multiple possible remappings
    mult_map_counts = pd.Series(
        pd.unique(interaction_dataframe[columns].values.ravel(
            'K'))).str.contains('MULT_MAP').sum()
    mult_map_interaction_count = interaction_dataframe.loc[interaction_dataframe[columns[0]].str.contains('MULT_MAP') | interaction_dataframe[columns[1]].str.contains('MULT_MAP')].shape[0]

    print(
        '{} identifiers were be remapped to multiple UniProt accession numbers (denoted by the prefix "MULT_MAP" in the dataset) for a total of {} interactions.\n'.
        format(mult_map_counts, mult_map_interaction_count))

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
    These identifiers should be removed because otherwise this would
    artificially inflate their numbers during frequent item set mining.

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

    print(
        'Omitted {} PPIs due to the existance of multiple mappings.\n'.format(
            np.sum(~selection)))
    interaction_dataframe = interaction_dataframe.loc[selection].apply(
        lambda x: x)  # TODO: fix this?
    interaction_dataframe = interaction_dataframe.reset_index(drop=True)
    return interaction_dataframe


# def lookup(interaction_dataframe, columns=None):
#     """ Remap identifiers to UniProt ACs.
#
#     Tries to replace non-UniProt ACs for the specified columns of a pandas DataFrame by looking in the unique
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


def check_unique_identifier(interaction_dataframe, columns=None):
    """
    Check, attempt to fix and report on whether the unique identifier column of
    a PSI-MITAB dataset is correct.

    E.g. HPDIB contains a few entries such as:
    dip:DIP-10051N|refseq:NP_285937
    where two identifiers are separated by a pipe symbol.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing psi-mitab protein-protein interactions.
    columns : list
        The names of the columns containing the unique identifiers.
        (The defaults are xref_A and xref_B).

    Returns
    -------
    None
        Modifies the supplied DataFrame in-place.
    """

    if not columns:
        columns = ['xref_A', 'xref_B']

    def lambda_helper(row_entry):
        split_entry = row_entry.split('|')
        uni_id = [x for x in split_entry if 'uniprot' in x]
        xref_id = [x for x in split_entry if 'refseq' in x]
        first_id = split_entry[0]
        if uni_id:
            if len(uni_id) != 1:
                print(
                    'WARNING: the interaction involving {} contains unresolvable interaction identifiers.\n'.
                    format(row_entry))
            new_id = uni_id[0]
        elif xref_id:
            if len(xref_id) != 1:
                print(
                    'WARNING: the interaction involving {} contains unresolvable interaction identifiers.\n'.
                    format(row_entry))
            new_id = xref_id[0]
        else:
            new_id = first_id
        return new_id

    for i in columns:
        non_unique_mask = interaction_dataframe[i].str.contains('\|')

        if any(non_unique_mask):
            interaction_dataframe.loc[non_unique_mask, i] = \
                interaction_dataframe.loc[non_unique_mask, i].apply(lambda_helper)
            print(
                '{} unique identifiers in column {} did not conform to the PSI-MITAB format (multiple identifiers were present) and a fix was attempted.\n'.
                format(len(non_unique_mask), i))

    # TODO: Find a way to monitor how many of the bad identifiers can be remapped to uniprot/refseq/other.
