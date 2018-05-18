#!/usr/bin/env python3
"""
Module to import protein-protein interaction networks in a uniform manner to a Pandas DataFrame.

The relevant columns are labelled as follows, as per the MI-TAB 2.5 format.

'xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B', 'aliases_A','aliases_B',
'detection_method', 'author', 'publication', 'taxid_A', 'taxid_B','interaction_type',
'source_database_ids', 'interaction_identifiers', 'confidence_score'

Intact:
1 #ID(s) interactor A
2 ID(s) interactor B
3 Alt. ID(s) interactor A
4 Alt. ID(s) interactor B
5 Alias(es) interactor A
6 Alias(es) interactor B
7 Interaction detection method(s)
8 Publication 1st author(s)
9 Publication Identifier(s)
10 Taxid interactor A
11 Taxid interactor B
12 Interaction type(s)
13 Source database(s)
14 Interaction identifier(s)
15 Confidence value(s)

hpidb:
1 # protein_xref_1
2 protein_xref_2
3 alternative_identifiers_1
4 alternative_identifiers_2
5 protein_alias_1
6 protein_alias_2
7 detection_method
8 author_name
9 pmid
10 protein_taxid_1
11 protein_taxid_2
12 interaction_type
13 source_database_id
14 database_identifier
15 confidence

BIOGRID:
1 #ID Interactor A
2 ID Interactor B
3 Alt IDs Interactor A
4 Alt IDs Interactor B
5 Aliases Interactor A
6 Aliases Interactor B
7 Interaction Detection Method
8 Publication 1st Author
9 Publication Identifiers
10 Taxid Interactor A
11 Taxid Interactor B
12 Interaction Types
13 Source Database
14 Interaction Identifiers
15 Confidence Values

virhost:
no header, same format

Additional columns are not omitted, but stored under their original header (if available).

Additional database specific notes:

- HPIDB2: Contains a header with 15/26 fields.
          The "_plus" dataset contains additional columns such as the protein sequence.
          Hosts and pathogens are stored in fixed columns, A and B respectively.
- VirHostNet2: No header available. Matches the number of columns.
               How many taxons?
               Hosts and pathogens are not stored in fixed columns (i.e. xref_A contains both hosts and pathogens).
               Column entries must be switched around to remedy this.
               Includes intra-species interactions. Mostly viral-viral?
               Identifiers are all UniProt ACs.
- PHISTO: Aberrant data format.
          Only contains pathogen-to-human interactions.
          Hosts and pathogens stored in fixed columns.
          No aliases.
          Identifiers are all UniProt ACs.

Some bash code to explore the datasets:

Number of columns in data set
$ head -1 hpidb2_March14_2017_mitab_plus.txt | awk -F'\t' '{print NF }'
26
What if there are some lines with more/fewer entries?
$ awk -F'\t' '{print NF}' hpidb2_March14_2017_mitab_plus.txt | sort -nu | tail -n 1
26
$ awk -F'\t' ' { for (i = 1; i <= NF; ++i) print i, $i; exit } ' hpidb2_March14_2017_mitab.txt
1 # protein_xref_1
2 protein_xref_2
3 alternative_identifiers_1
4 alternative_identifiers_2
5 protein_alias_1
6 protein_alias_2
7 detection_method
8 author_name
9 pmid
10 protein_taxid_1
11 protein_taxid_2
12 interaction_type
13 source_database_id
14 database_identifier
15 confidence

Retrieve data sources present in PPI dataset
$ tail -n +2 hpidb2_March14_2017_mitab.txt | cut -f1 | sed -r 's/(^.*):.*/\1/g' | sort -u
ensembl
entrez gene/locuslink
intact
uniprotkb

- check for overlap/redundancy
- check for inter-intra
- check for consistency in column A versus B
- check for number of hosts
- check for number of interactions per host-viral pair

"""

import pandas as pd

from pathlib import Path


def read_mi_tab(filepath, encoding='ISO-8859-1'):
    """Convert a protein-protein interaction file in PSI-MITAB 2.5 format
    into a pandas DataFrame.

    Parameters
    ----------
    filepath : str
        Filepath to mi-tab file.
    encoding : str
        Encoding of the file.

    Returns
    -------
    DataFrame
        pandas DataFrame containing the interaction data.
    """
    file = Path(filepath)
    with open(file, encoding=encoding) as f:
        first_line = f.readline()
    # Check if a header is present
    if any(i in first_line.lower() for i in [
            '#', 'Interactor A', 'protein_xref_1', 'author_name',
            'source database'
    ]):
        header = 0
    else:
        header = None
    # read file into pandas, setting encoding and header
    df = pd.read_csv(file, sep='\t', encoding=encoding, header=header)
    # because there is no standard set of column names across mitab files from
    # different sources, change them to following
    col_names = [
        'xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B',
        'aliases_A', 'aliases_B', 'detection_method', 'author', 'publication',
        'taxid_A', 'taxid_B', 'interaction_type', 'source_database_ids',
        'interaction_identifiers', 'confidence_score'
    ]
    # if a file contains more than the standard 15 columns, e.g. hpidb-plus
    # only rename the first 15 ones and keep the rest
    df.columns = col_names + df.columns.tolist()[len(col_names):]
    # homogenise data entry formats e.g. for the hpidb2 format style
    # taxid:9606(human|Homo sapiens)
    # to
    # taxid:9606
    df.taxid_A = df.taxid_A.str.split('(').str.get(0)
    df.taxid_B = df.taxid_B.str.split('(').str[0]
    # name dataframe
    df.name = file.name
    print('Read PPI data set', df.name, 'from', str(filepath))

    return df


def read_psi_mi_tab(filepath, name):
    """Read HPIDB2 data set into pandas DataFrame.

    Renames headers to standard format used by this module.
    Can deal with both the standard and -plus data set provided by the HPIDB2.

    Parameters
    ----------
    filepath : string or Pathlib object
        file path to HPIDB2 data set.

    Returns
    -------
    DataFrame
        DataFrame containing the HPIDB2 interaction data.

    """
    # TODO: notification for regular or plus file...
    file = Path(filepath)
    df = pd.read_csv(file, sep='\t', encoding='ISO-8859-1', header=0)
    new_partial_columns = [
        'xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B',
        'aliases_A', 'aliases_B', 'detection_method', 'author', 'publication',
        'taxid_A', 'taxid_B', 'interaction_type', 'source_database_ids',
        'interaction_identifiers', 'confidence_score'
    ]
    # re-name first few columns to same format used for VirHost data set.
    df.columns = new_partial_columns + df.columns.tolist()[len(
        new_partial_columns):]
    # rename_dict = dict(zip(df.columns.tolist()[:len(new_partial_columns)], new_partial_columns))
    # df.rename(columns=rename_dict, inplace=True)
    df.taxid_A = df.taxid_A.str.split('(').str.get(0)
    df.taxid_B = df.taxid_B.str.split('(').str[0]
    # name dataframe
    df.name = name
    print('Read PPI data set', df.name, 'from', str(filepath) + '.')
    return df


def read_mitab_virhost(filepath):
    """Read VirHost data set into pandas DataFrame.

    Renames headers to standard format used by this module.

    Parameters
    ----------
    filepath : string or Pathlib object
        file path to VirHost data set.

    Returns
    -------
    DataFrame
        DataFrame containing the VirHost interaction data.

    """
    file = Path(filepath)
    df = pd.read_csv(
        file,
        sep='\t',
        encoding='ISO-8859-1',
        names=[
            'xref_A', 'xref_B', 'alt_identifiers_A', 'alt_identifiers_B',
            'aliases_A', 'aliases_B', 'detection_method', 'author',
            'publication', 'taxid_A', 'taxid_B', 'interaction_type',
            'source_database_ids', 'interaction_identifiers',
            'confidence_score'
        ])
    df.detection_method = df.detection_method.str.replace('"', '')
    df.interaction_type = df.interaction_type.str.replace('"', '')
    df.source_database_ids = df.source_database_ids.str.replace('"', '')
    # name dataframe
    df.name = 'VirHostNet2'
    print('Read PPI data set', df.name, 'from', str(filepath) + '.')
    return df


def read_mitab_phisto(mitab_filepath, mi_filepath):
    """Read PHISTO data set into pandas DataFrame.

    Renames headers to standard format used by this module and modifies a
    number of entries to match the standard format.

    Parameters
    ----------
    mitab_filepath : string
        file path to PHISTO data set.
    mi_filepath : string or Pathlib object
        file path to psi-mi .obo file.

    Returns
    -------
    DataFrame
        DataFrame containing the PHISTO interaction data.

    """
    file = Path(mitab_filepath)
    # display_id = entry name on uniprot
    # entry B is always pathogen partner in PHISTO
    df = pd.read_csv(
        file,
        sep=',',
        encoding='ISO-8859-1',
        header=0,
        names=[
            'pathogen', 'taxid_B', 'xref_B', 'display_id_B', 'xref_A',
            'display_id_A', 'detection_method', 'publication'
        ])
    # Retrieve molecular interaction obo ontology
    mi_dict = phisto_load_mi_ontology(mi_filepath)
    reverse_mi_dict = {name: id for id, name in mi_dict.items()}
    # Fetch psi-mi id for given detection method in PHISTO interaction file,
    # otherwise, keep the original detection method name. Required because of e.g. "Other Methods".
    df.detection_method = df.detection_method.map(
        lambda x: 'psi-mi:' + str(reverse_mi_dict.get(x, x)) + '(' + str(x) + ')'
    )
    # Append column name in front of entries
    df.publication = df.publication.map(lambda x: 'pubmed:' + str(x))
    df.xref_A = df.xref_A.map(lambda x: 'uniprotkb:' + str(x))
    df.xref_B = df.xref_B.map(lambda x: 'uniprotkb:' + str(x))
    df.taxid_B = df.taxid_B.map(lambda x: 'taxid:' + str(x))
    # Add Human taxid_A column
    df['taxid_A'] = 'taxid:9606'
    # name dataframe
    df.name = file.name
    print('Read PPI data set', df.name, 'from', str(mitab_filepath))
    return df


def phisto_load_mi_ontology(filepath):
    """ Read in PSI-MI molecular interaction ontology file and store in dictionary.

    Ontology file should be in .obo file.
    Can be retrieved from http://ontologies.berkeleybop.org/mi.obo


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
