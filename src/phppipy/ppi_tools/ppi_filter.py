import numpy as np
import pandas as pd


def annotate_inter_intra_species(interaction_dataframe):
    """Adds column to DataFrame specifying whether interaction is
    inter- or intra-species.

    The added column is named "inter-intra-species".

    Parameters
    ----------
    interaction_dataframe : DataFrame
        The pandas DataFrame should correspond to the PSI-MITAB format.

    Returns
    -------
    None
        Modifies DataFrame in-place by adding the "inter-intra-species" column.

    """
    interaction_dataframe['inter-intra-species'] = np.where(
        interaction_dataframe.taxid_A == interaction_dataframe.taxid_B,
        'intra', 'inter')


def annotate_inter_intra_pathogen(interaction_dataframe, pathogen_taxa):
    """Adds column to DataFrame specifying whether interaction is between
    pathogens or between pathogens and hosts.

    The added column is named "inter-intra-pathogen".

    Parameters
    ----------
    interaction_dataframe : DataFrame
        The pandas DataFrame should correspond to the PSI-MITAB format.
    pathogen_taxa : array-like
        A list of pathogen taxon IDs of the form: ['taxid:10377', ...]

    Returns
    -------
    None
        Modifies DataFrame in-place by adding the "inter-intra-pathogen" column.

    """
    interaction_dataframe['inter-intra-pathogen'] = np.where(
        interaction_dataframe.taxid_A.isin(pathogen_taxa) &
        interaction_dataframe.taxid_B.isin(pathogen_taxa), 'intra', 'inter')


def unique_identifier(df, column_name='xref_partners_sorted', columns=None):
    """ Creates a new column containing a sorted string of two
    identifiers found in the supplied columns, i.e. a unique
    reference identifier for an association pair.

    Parameters
    ----------
    df : DataFrame
        Protein-protein interaction dataframe.
    column_name : str
        The name of the new column.
    columns : list
        The columns to use for the creation of a new unique identifier.
        Must contain exactly two columns.

    Returns
    -------
    None
        Modifies DataFrame in-place by adding the unique identifier column.
    """

    if not columns:
        columns = ['xref_A', 'xref_B']

    xref_partners_sorted_array = np.sort(
        np.stack((df.xref_A, df.xref_B), axis=1), axis=1)

    xref_partners_df = pd.DataFrame(
        xref_partners_sorted_array, columns=['A', 'B'])

    df[column_name] = xref_partners_df['A'] + '%' + xref_partners_df['B']


def reorder_pathogen_host_entries(interaction_dataframe, host_list):
    """ Moves all pathogen entries to B columns and host entries to A columns.

    Selects all interaction entries where host entries occur in B columns instead of A and swaps the A and B columns.
    Or vice versa, where pathogen occurs in column A.

    # Note: HPIDB2 always has hosts as partner A and PHISTO has human/pathogen labeled as well.

    https://stackoverflow.com/questions/25792619/what-is-correct-syntax-to-swap-column-values-for-selected-rows-in-a-pandas-data

    Parameters
    ----------
    interaction_dataframe : DataFrame
        The pandas DataFrame containing the protein-protein interactions that need to be sorted.
    host_list : list
        List of host taxonids to look for in the B column.

    Returns
    -------
    None
        Modifies the interaction DataFrame inplace.
    """

    host_list = set(host_list)
    host_position_mask = interaction_dataframe['taxid_B'].isin(host_list)
    column_names = [
        'xref', 'taxid', 'aliases', 'alt_identifiers', 'display_id',
        'taxid_name'
    ]
    columns_to_swap = [
        name + label for name in column_names for label in ['_A', '_B']
    ]
    columns_after_swap = [
        name + label for name in column_names for label in ['_B', '_A']
    ]
    interaction_dataframe.loc[
        host_position_mask, columns_to_swap] = interaction_dataframe.loc[
            host_position_mask, columns_after_swap].values
