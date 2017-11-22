#!/usr/bin/env python3

"""
Main module to parse protein-protein interaction datasets (mitab format), label them with GO and InterPro annotations
and convert them into transaction data sets suitable for frequent item set mining.

Must be run as a script.

Paths are hard-coded relative to the repository.
"""

# TODO: mutuable list default values should be changed to tuples or if None, set to list

import os, sys
sys.path.append(os.path.abspath('..'))
# or Path('.').parent

from pathlib import Path

import argparse
import numpy as np
import pandas as pd

import id_mapper
import label_go
import label_interpro
import ppi_import

from data_prep import retrieve_taxids

from go_tools import gaf_parser
from go_tools import obo_tools

# import statements for interactive console
# import sys
# sys.path.insert(0, '/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/go-tools/go-tools')
# from go_tools import gaf_parser
# from go_tools import obo_tools
# sys.path.append("/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi_scripts/go_enrichment_tool")
# sys.path.append("/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi_scripts")

# pandas display options
pd.set_option('display.max_columns', None)
pd.options.display.max_colwidth = 200


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

    concatenated_df = pd.concat(list_of_datasets, axis=0, ignore_index=True, join='outer')
    concatenated_df.reset_index(inplace=True, drop=True)
    # df_concat = df_concat.drop_duplicates(subset=['xref_partners_sorted'])

    return concatenated_df


def reorder_pathogen_host_entries(interaction_dataframe, host_list=list('taxid:9606')):
    """ Moves all pathogen entries to B columns and host entries to A columns.

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
    """
    host_position_mask = interaction_dataframe['taxid_B'].isin(host_list)
    # Note: HPIDB2 always has hosts as partner A and PHISTO has human/pathogen labeled too.
    column_names = ['xref', 'taxid', 'aliases', 'alt_identifiers']
    columns_to_swap = [name + label for name in column_names for label in ['_A', '_B']]
    columns_after_swap = [name + label for name in column_names for label in ['_B', '_A']]
    interaction_dataframe.loc[host_position_mask, columns_to_swap] = interaction_dataframe.loc[host_position_mask,
                                                                                               columns_after_swap].values


def label_host_pathogen(interaction_dataframe, pathogen_set, columns=list(('xref_A_GO', 'xref_B_GO')),
                        taxid_columns=list(('taxid_A', 'taxid_B'))):
    """ Adds virus/host label to selected columns in DataFrame.

    The columns and taxid lists should be of equal length and ordered in the same manner. E.g. for each position
    in the columns list, the corresponding position in the taxid list should repeat the associated taxid column.
    Calls the label_term() function which takes a set as its input and returns a string with comma separated terms.

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
    """
    if not len(columns) == len(taxid_columns):
        raise ValueError("Lists of columns to label should match taxid columns.")

    else:
        for col, taxid in zip(columns, taxid_columns):
            lambda_label = lambda x: label_term(x[col], x[taxid], pathogen_set)
            interaction_dataframe[col] = interaction_dataframe.apply(lambda_label, axis=1)


def label_term(label_set, taxid, pathogen_set):
    """Adds host/virus labels to terms and returns string representation (comma separated terms).

    Should be called via an apply function that calls it for both partners' term sets (e.g. GO terms) and taxids
    for each interaction in a dataset.

    Note: currently all hosts are labelled identically because the function simply checks whether a taxid belongs
          to a set of virus taxids.


    Parameters
    ----------
    label_set : set
        The set of terms associated with a specific protein partner, e.g. those stored in the xref_A/B_GO
        columns after running the annotate_GO() function on an interaction DataFrame.
    taxid : string
        The taxid of the protein partners, stored in the taxid_A/B columns of the interaction DataFrame.
        E.g. taxid:9606
    pathogen_set : set
        A set containing all relevant pathogen taxids.

    Returns
    -------
    string
        A string of labeled terms separated by commas.

    """
    if pd.isnull(label_set):
        return label_set
    else:
        if taxid in pathogen_set:
            label = 'v@'
        else:
            label = 'h@'

        labeled_list = []
        for i in label_set:
            labeled_list.append(label + i)
        return ','.join(labeled_list)


def pathogen_group_mapper(taxid, pathogen_group_dict):
    """Map pathogen taxid to higher order taxonomic grouping.

    Parameters
    ----------
    taxid : str
        A pathogen taxid of the format "taxid:#####".
    pathogen_group_dict : dictionary
        A dictionary mapping pathogen groups to taxids.
        E.g. pathogen_group_dict['bovine_ah1'] =
        ['10320', '10322', '79889', '79890', '10323', '10324', '31517', '31518', '10321', '31519', '45407',
        '1548190', '1548191', '1548192']

    Returns
    -------
    str
        The pathogen taxonomic group name.

    """
    for pathogen_group, taxids in pathogen_group_dict.items():
        if taxid in taxids:
            return pathogen_group
    return np.NaN

    # lambda_pathogen_mapper = lambda x: [pathogen_group for pathogen_group, taxids in pathogen_group_dict.items()
    #                                     if taxid in taxids]


if __name__ == '__main__':
    # Check provided arguments
    parser = argparse.ArgumentParser(
        description='Script to filter and annotate PPI databases into a dataset suitable for itemset mining.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dry', action="store_true", help='Perform dry run without saving output.')
    parser.add_argument('-o', '--output', type=str, default=r'../../data/processed/transaction_datasets/',
                        help='Output file or path')
    # TODO: add arguments for every type of input data source instead of hardcoding directory structure.
    args = parser.parse_args()

    # Read in PPI datasets
    df_virhost = ppi_import.read_mitab_virhost(r'../../data/raw/ppi_data/VirHostNet_January_2017.txt')

    df_hpidb2 = ppi_import.read_mitab_hpidb2(r'../../data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt')

    df_phisto = ppi_import.read_mitab_phisto(r'../../data/raw/ppi_data/phisto_Jan19_2017.csv',
                                             r'../../data/raw/ppi_data/mi.obo')

    # Concatenate the different sources
    print('Concatenating PPI datasets...')
    df_concat = concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto])

    # Map non-UniProt entries to UniProt ACs
    id_mapper.map2uniprot(df_concat, filepath=r'../../data/interim/mappings/')
    # TODO: no intact-EBI mapping...
    # Remove PPIs with multiple id mappings because otherwise this would artificially inflate their numbers during
    # frequent item set mining
    df_concat = id_mapper.remove_mult(df_concat)

    # Add unique identifier for interaction pairs
    # note: absolutely critical to have index from 0 to n here!
    xref_partners_sorted_array = np.sort(np.stack((df_concat.xref_A, df_concat.xref_B), axis=1), axis=1)
    # df_concat['xref_partners_sorted'] = pd.Series(map(tuple, xref_partners_sorted_array))
    xref_partners_df = pd.DataFrame(xref_partners_sorted_array, columns=['A', 'B'])
    df_concat['xref_partners_sorted'] = xref_partners_df['A'] + '%' + xref_partners_df['B']

    # Label interactions as being within or between species.
    annotate_inter_intra(df_concat)

    # Filter out intra-species interactions
    print('Omitting intra-species interactions...')
    df_concat = df_concat[df_concat['inter-intra'] == 'inter-species']
    df_concat = df_concat.reset_index(drop=True)

    # Size of each data source
    print('\nData source sizes:\n')
    print(df_concat.groupby('origin').size())

    # Count duplicates before any filtering
    # TODO: check if duplicates from one data source differ in detection method, publication or something else
    print('\nNumber of duplicated interactions on raw datasets')
    print(np.sum(df_concat.duplicated(subset=['xref_partners_sorted'])))

    print('\nNumber of unique interactions per raw data set')
    print(df_concat.groupby('origin')['xref_partners_sorted'].nunique())
    # bash script to check overlap:
    # comm <(cut -f3 -d, phisto_Jan19_2017.csv | sed 's/"//g' | sort -u ) <(cut -f2 hpidb2_March14_2017_mitab.txt | sed s/uniprotkb://g | sort -u)

    # Check non-UniProt ACs
    print('\nNumber of interactions without UniProt AC')
    print(df_concat.loc[(~df_concat.xref_A.str.contains('uniprot', case=False)) |
                        (~df_concat.xref_B.str.contains('uniprot', case=False))].groupby('origin').size())

    # Remove duplicate interaction pairs (including different detection methods and publications)
    # https://stackoverflow.com/a/41650846
    # https://stackoverflow.com/questions/33042777/
    # Note that this will result in the first dataset (e.g. hpidb2) having priority over the others.
    df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted'], keep='first')
    # df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted', 'taxid_B', 'taxid_A'], keep='first')
    df_concat_dedup = df_concat_dedup.reset_index(drop=True)
    # When subsetting duplicates, also check taxid_A and taxid_B
    # e.g. "Human herpesvirus 1 STRAIN KOS","10306","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
    # "Human herpesvirus 1","10298","A1Z0Q5","GD_HHV1K ","Q15223","NECT1_HUMAN ","fluorescence-activated cell sorting","12011057"
    # TODO: this could result in duplicates for item set mining...
    # These differ in pathogen strain, but it's the same interaction.
    # Will come into play when comparing groups and different groups have same interaction
    # Which is not the case here, since both are human herpesvirus 1 strains
    # Similar issue for pubmed ID, method, etc.
    print('\nOmitting duplicates as defined by UniProt ACs and taxids.')

    # Retrieve Herpesviridae (taxid:10292) list, see retrieve_taxids.py script to generate child taxids
    # TODO: import from retrieve_Taxids and create on the spot
    # TODO: then combine this code into a function to subset a given dataframe for a given taxid and its children
    # TODO: and a host list or default to all hosts
    # TODO: note that this filtering will generally result in intra-host interactions being omitted,
    # TODO: while retaining intra-viral ones
    with Path(r'../../data/interim/child_taxids_of_10292.txt').open() as taxid_file:
        herpes_taxids = [str('taxid:' + line.split('|')[0]) for line in taxid_file]
    print(r'Retrieving Herpes taxids from ../../data/interim/child_taxids_of_10292.txt')

    # Extract Herpesviridae interactions
    # This omits all human-human interactions, but not intra-virus interactions
    df_herpes = df_concat_dedup.loc[(df_concat_dedup.taxid_A.isin(herpes_taxids)) |
                                    df_concat_dedup.taxid_B.isin(herpes_taxids)]
    df_herpes = df_herpes.reset_index(drop=True)
    print('Omitting all non-Herpes interactions.')

    # Check how many non-UniProt interactions are left
    print('\nNumber of non-UniProt AC interactions for Herpes interactions')
    print(df_herpes.loc[(~df_herpes.xref_A.str.contains('uniprot')) |
                        (~df_herpes.xref_B.str.contains('uniprot'))].groupby('origin').size())

    # Filter on UniProt ACs due to GO mapping being UniProt-based
    df_herpes = df_herpes.loc[(df_herpes.xref_A.str.contains('uniprot')) & (df_herpes.xref_B.str.contains('uniprot'))]
    df_herpes = df_herpes.reset_index(drop=True)
    print('Omitting all non-UniProt AC entries...')

    # TODO: Filter on interaction types
    # psi - mi: MI:0915(physical association)
    # psi - mi: MI:0914(association)
    # psi - mi: MI:0407(directinteraction)
    # psi - mi: MI:0403(colocalization)
    # psi - mi: MI:0217(phosphorylation reaction)
    # psi - mi: MI:0208(geneticinteraction)
    # print(df_herpes['interaction_type'].value_counts())
    # df_herpes = df_herpes.loc[df_herpes.interaction_type.str.contains('association')] # fails due to NaNs... df.isnull().sum()


    # Print some information about the interactions
    print('\nNumber of inter versus intra interactions:')
    print(df_herpes.groupby(['inter-intra', 'origin']).size())
    print('\nNumber of taxid involved in intra-species interactions')
    print(df_herpes.loc[df_herpes['inter-intra'] == 'intra-species'].groupby(['taxid_A', 'taxid_B']).size())  # no human
    print('\nNumber of human-human interactions')
    print(df_herpes.loc[(df_herpes['inter-intra'] == 'intra-species') &
                        (df_herpes['taxid_A'].str.contains('9606'))].shape)
    print('\nNumber of interactions without UniProt AC')
    print(df_herpes.loc[(~df_herpes.xref_A.str.contains('uniprot')) |
                        (~df_herpes.xref_B.str.contains('uniprot'))].shape)

    # Check which hosts are present
    print('\nNumber of interactions per host')
    all_taxids = df_herpes['taxid_A'].append(df_herpes['taxid_B']).unique()
    host_taxids = list(np.setdiff1d(all_taxids, herpes_taxids))

    taxid_names_path = Path(r'../../data/raw/taxdump/names.dmp')
    name2taxid, taxid2name = retrieve_taxids.parse_taxid_names(str(taxid_names_path))

    for i in host_taxids:
        taxid = i.split(':')[1]
        count = df_herpes['xref_partners_sorted'].loc[(df_herpes['taxid_A'] == i) | (df_herpes['taxid_B'] == i)].shape
        print(taxid, taxid2name[taxid], count)

    # Create taxid dictionaries
    taxid_nodes_path = Path(r'../../data/raw/taxdump/nodes.dmp')
    taxid2parent, taxid2rank = retrieve_taxids.parse_taxid_nodes(str(taxid_nodes_path))
    parent2child = retrieve_taxids.create_parent2child_dict(taxid2parent)

    # Move all host partners in xref_B to xref_A (only an issue for VirHostNet)
    # Note, also move ALL other associated labels...
    reorder_pathogen_host_entries(df_herpes, host_taxids)

    # Count missing values across columns
    print('\nMissing values in each column:')
    print(df_herpes.isnull().sum(axis=0), '\n')

    # Create Gene Ontology dictionaries
    print('Creating Gene Ontology dictionaries...')
    # TODO: create quickGO query using host list + herpesviridae
    go_dict = obo_tools.importOBO(r'../../data/raw/go_data/go.obo')
    # TODO: allow selection of GO namespaces
    obo_tools.buildGOtree(go_dict, root_nodes=['GO:0008150', 'GO:0005575', 'GO:0003674'])

    # protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.split(':').str[1].unique())
    # Needs special extract pattern because of PTM processing id's e.g. P04295-PRO_0000115940
    protein_set = set(df_herpes.xref_A.append(df_herpes.xref_B, ignore_index=True).str.extract('^.*:(\w*)-?',
                                                                                               expand=False).unique())
    gaf_dict = gaf_parser.importGAF(r'../../data/raw/go_data/gene_association_hosts_10292.goa', protein_set)

    # Number of proteins without GO annotation
    print('\nNumber of proteins without GO annotation')
    not_annotated = [i for i in protein_set if i not in gaf_dict]
    print(len(not_annotated))

    # Add GO annotations
    label_go.annotate_GO(df_herpes, gaf_dict)

    # Re-map GO depth
    # TODO: if part_of relations are allowed, nodes are connected across ontologies.
    # TODO: either label ontologies or keep them all separated?
    # TODO: requires a check during mapping to see if terms have the same namespace or not!
    depth = {'biological_process': 2, 'molecular_function': 1, 'cellular_component': 1}
    exclusion = {'GO:0005623': 'cell', 'GO:0005488': 'binding', 'GO:0043226': 'organelle',
                 'GO:0044422': 'organelle part', 'GO:0044464': 'cell part', 'GO:0033643': 'host cell part',
                 'GO:0033646': 'host intracellular part', 'GO:0043656': 'intracellular region of host',
                 'GO:0043657': 'host cell', 'GO:0018995': 'host', 'GO:0044424': 'intracellular part',
                 'GO:0016032': 'viral process', 'GO:0044215': 'other organism',
                 'GO:0050789': 'regulation of biological process', 'GO:0005515': 'protein binding',
                 'GO:0019012': 'virion', 'GO:0044423': 'virion part', 'GO:0039642': 'virion nucleoid',
                 'GO:0019028': 'viral capsid', 'GO:0055036': 'virion membrane', 'GO:0036338': 'viral membrane',
                 'GO:0098015': 'virus tail'}
    # np.sum(df_herpes1.loc[df_herpes1['xref_A_GO'].notnull(), 'xref_A_GO'].apply(lambda x: True if 'GO:0043657' in x else False))
    # regulation of biological process: 5108 and 2929 in A/B
    # cell part 6827 and 1
    # organelle part 5317 and 1
    # organelle 6203 and 0
    # binding 6607 and 3224
    # protein binding 5518 and 2442
    # np.sum(df_herpes.loc[df_herpes['xref_B_GO'].notnull(), 'xref_B_GO'].apply(lambda x: True if any([i in x for i in go_dict.get('GO:0005515').recursive_children]) else False))
    # only 186 of protein binding child terms annotated for virus (4100 for host)
    # cell 23 and 0
    # virion 91 and 1509
    # none of child terms present..
    # virion part 92 and 1456
    # 1456 virus child terms, 92 host

    print('Remapping GO terms to desired depth', depth)
    label_go.remap_GO_depth(df_herpes, depth, go_dict, list(exclusion))


    # depth = {'biological_process': 2, 'molecular_function': 1, 'cellular_component': 2}
    # remap_GO_depth(df_herpes, depth, go_dict, ['xref_B_GO'])
    # MF: needs deeper mapping for "binding"
    # CC: OK, virus needs deeper...membrane would benefit from deeper as well.

    def count_term(go_term, interaction_dataframe, column_list=None):
        if column_list is None:
            column_list = ['xref_A_GO', 'xref_B_GO']
        count = 0
        sets_array = interaction_dataframe[column_list].values.ravel()
        sets_array = sets_array[~pd.isnull(sets_array)]
        for i in sets_array:
            if go_term in i:
                count += 1
        print(f'Count for {go_dict.get(go_term).name}')
        return count


    # Add virus/host labels to GO and converts set to string of comma separated values
    label_host_pathogen(df_herpes, herpes_taxids, columns=['xref_A_GO', 'xref_B_GO'],
                        taxid_columns=['taxid_A', 'taxid_B'])

    # Map pathogen id to higher label
    pathogen_group_dict = {'bovine_ah1': '10320', 'bovine_hv1': '79889', 'epstein_barr': '10376', 'equid_av1': '10326',
                           'equid_gv2': '12657', 'gallid_av2': '10390', 'human_hsv1': '10298', 'saimiri_gv2': '10381',
                           'human_av2': '10310', 'human_av3': '10335', 'human_bv5': '10359', 'human_gv8': '37296',
                           'human_bv6A': '32603', 'human_bv6B': '32604', 'murid_bv1': '10366', 'murid_gv4': '33708',
                           'papiine_gv1': '106332', 'suid_av1': '10345'}

    for i, j in pathogen_group_dict.items():
        pathogen_group_dict[i] = [j] + retrieve_taxids.get_children(j, parent2child)

    df_herpes['pathogen_groups'] = df_herpes.apply(lambda x: pathogen_group_mapper(x['taxid_B'].split(':')[1],
                                                                                   pathogen_group_dict), axis=1)

    print('\nNumber of interactions for each pathogen grouping')
    print(df_herpes.groupby('pathogen_groups').size())

    # unique_ac = set(df_herpes.xref_A.append(df_herpes.xref_B,
    # ignore_index=True).str.extract('^.*:(\w*)-?',expand=False).unique())

    # Get set of unique ACs
    unique_ac = set(pd.unique(df_concat['xref_B'].str.extract('^.*:(\w*)-?', expand=False).append(
        df_concat['xref_A'].str.extract('^.*:(\w*)-?', expand=False), ignore_index=True)))

    # Add InterPro labels
    print('Adding InterPro annotations...')
    uniprot2interpro = label_interpro.create_uniprot2interpro_dict(unique_ac,
                                                                   filepath=r'../../data/interim/interpro_data'
                                                                            r'/protein2ipr_filtered.txt')
    label_interpro.annotate_interpro(df_herpes, uniprot2interpro)
    label_host_pathogen(df_herpes, herpes_taxids, columns=['interpro_A', 'interpro_B'],
                        taxid_columns=['taxid_A', 'taxid_B'])

    # Save to transaction database
    # WRONG APPROACH: NaN values in either column will result in NaN in combined column...
    # combined_GO_labels_old = df_herpes['xref_A_GO'] + ',' + df_herpes['xref_B_GO']
    # compare combined_GO_labels[[8061, 8062, 8122,8123]] with combined_GO_labels_old[[8061, 8062, 8122,8123]]
    # and df_herpes.loc[[8061, 8062, 8122,8123],['xref_A_GO','xref_B_GO']]
    # see https://stackoverflow.com/questions/26614465/python-pandas-apply-function-if-a-column-value-is-not-null

    # Merge A and B GO columns into one
    has_GO_mask = df_herpes[['xref_A_GO', 'xref_B_GO']].notnull().all(axis=1)
    combined_GO_labels = df_herpes.loc[has_GO_mask, 'xref_A_GO'] + ',' + df_herpes.loc[has_GO_mask, 'xref_B_GO']
    # updates NaN in called Series/DF with values from argument Series/DF
    # also supplies values for non-existing indices/columns
    # e.g. indices with NaNs were absent from combined_GO_labels Series and are now added again.
    combined_GO_labels = combined_GO_labels.combine_first(df_herpes['xref_A_GO'])
    combined_GO_labels = combined_GO_labels.combine_first(df_herpes['xref_B_GO'])
    # name Series and join into dataframe
    combined_GO_labels.rename('GO', inplace=True)
    df_herpes = df_herpes.join(combined_GO_labels)

    # Merge A and B InterPro columns into one
    has_InterPro_mask = df_herpes[['interpro_A', 'interpro_B']].notnull().all(axis=1)
    combined_interpro_labels = df_herpes.loc[has_InterPro_mask, 'interpro_A'] + ',' + df_herpes.loc[
        has_InterPro_mask, 'interpro_B']
    # updates NaN in called Series/DF with values from argument Series/DF
    combined_interpro_labels = combined_interpro_labels.combine_first(df_herpes['interpro_A'])
    combined_interpro_labels = combined_interpro_labels.combine_first(df_herpes['interpro_B'])
    # name Series and join into dataframe
    combined_interpro_labels.rename('interpro', inplace=True)
    df_herpes = df_herpes.join(combined_interpro_labels)

    # def combine_terms_for_export(interaction_dataframe, column_list=list(('xref_A_GO', 'xref_B_GO')), sep=','):
    #     # Combining more than two labels?
    #     # Difficulties: for loop, but leading and trailing comma? Use .str.strip afterwards?
    #     # df.combine_first() called for all combinations?
    #     # Create new column initialised to empty strings
    #     interaction_dataframe['combined_labels'] = ''
    #     # combined_labels = pd.DataFrame({'combined_labels':['' for i in range(interaction_dataframe[column_list].shape[0])]})
    #     # Append each term column to the combination column
    #     for i in column_list:
    #         interaction_dataframe['combined_labels'] += interaction_dataframe[i] + sep
    #     # NaN's in ANY column will reset the combined column to NaN again...
    #     #TODO: map()? or keep as list as long as possible and only then convert to string
    #     #TODO: label_set should then return a set object, and another function should be created to convert this to a string
    #     for i in column_list:
    #         combined_labels = combined_labels.combine_first(interaction_dataframe[i])
    #     combined_labels.rename('combined_labels', inplace=True)
    #     interaction_dataframe = interaction_dataframe.join(combined_labels)
    #
    #
    # combine_terms_for_export(df_herpes)

    # Save output
    if not args.dry:
        # Note: Pathlib functionality is broken in Pandas 0.20!
        output_directory = Path(args.output)
        output_directory.mkdir(exist_ok=True)
        output_path = output_directory / 'ppi_transactions.csv'
        print('Saving labelled PPI datasets to', output_directory.resolve())
        df_output = df_herpes[['xref_partners_sorted', 'GO', 'interpro']]
        df_output.reset_index(inplace=True, drop=True)
        df_output.to_csv(output_path, sep=',', index=False, header=False)
        # still requires script to remove quotes...
        with Path(output_path).open("r") as source:
            data = source.read()
            no_quotes = data.replace('"', '')
        with Path(output_path).open("w") as target:
            target.write(no_quotes)


        # optional: expand each label into its own column, but introduces a lot of empty columns
        # df_output = pd.concat([df_output['xref_partners_sorted'],
        #                        df_output['GO'].str.split(';', expand=True),
        #                        df_output['interpro'].str.split(';', expand=True)], axis=1)
        # df_output.to_csv(str(output_directory) + r'/ppi_transactions.csv', sep=',', index=False, header=False, doublequote=False)

        # Save to separate transaction datasets for each pathogen group
        # for i in df_herpes['pathogen_groups'].dropna().unique():
        for i in pd.unique(df_herpes['pathogen_groups'].dropna()):
            df_output_grouped = df_herpes.loc[df_herpes['pathogen_groups'] == i,
                                              ['xref_partners_sorted', 'GO', 'interpro']]
            df_output_grouped.reset_index(inplace=True, drop=True)
            df_output_grouped.to_csv(str(output_directory) + r'/' + str(i) + '.csv', sep=';', index=False, header=False)

        # TODO: df_herpes[['xref_A_GO', 'xref_B_GO']].notnull() how to melt this to 1 column of boolean indices?
        # df_herpes[['xref_A_GO', 'xref_B_GO']].notnull().all(axis=1)
        # TODO: EBI-intact identifiers?

        df_herpes.to_csv(str(output_directory) + r'/ppi_network.csv', sep=';', index=False)

# TODO: create taxid-pair identifier
