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
        interaction_dataframe.taxid_B.isin(pathogen_taxa), 'intra',
        'inter')
