#!/usr/bin/env python3
"""
Module to find pairwise associations between the annotation terms belonging
to protein interactions.
"""

import itertools
import math
import sys

import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.stats.multitest

from goscripts import obo_tools
from goscripts import gaf_parser

from joblib import Parallel, delayed


def merge_annotations(interaction_dataframe, columns):
    """Merges the annotations for a protein-protein interaction.

    Given two columns of a protein-protein interaction dataframe, each of
    which contains a type of annotation data, this function returns the merged
    set of those annotations.

    E.g. Column 1: {IPR011333, IPR003131, IPR000210}
         Column 2: {GO:0046872}
         Result: {GO:0046872, IPR011333, IPR003131, IPR000210}

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein-protein interactions and annotations.
    columns : list
        A list of annotation columns to merge. Expects the column contents to
        be array-like, not strings.

    Returns
    -------
    pandas Series of sets
        Array containing a set of annotations for each interaction (row).
    """

    # replace NaNs with empty sets
    for i in columns:
        interaction_dataframe.loc[interaction_dataframe[
            i].isnull(), i] = interaction_dataframe.loc[interaction_dataframe[
                i].isnull(), i].apply(lambda x: set())

    # join all sets in each supplied column
    # the unpacking operator can accept an array of lists or sets
    merged_annotations = interaction_dataframe[columns].apply(
        lambda x: set().union(*x), axis=1)
    """ Example usage
    Row entry:
        interpro_xref_A               {IPR029199, IPR026667}
        interpro_xref_B    {IPR027417, IPR001889, IPR013672}
        Name: 1, dtype: object
    Lambda result:
        {'IPR001889', 'IPR013672', 'IPR026667', 'IPR027417', 'IPR029199'}
    """

    return merged_annotations


def add_hp_label(merged_annotations_column, label_type):
    """Adds prefix to annotation labels that identify the annotation as
    belonging to the provided label_type (e.g. 'h@' for host proteins).

    Parameters
    ----------
    merged_annotations_column : array-like (pandas Series))
        An array containing sets of annotations that need to be labeled.
        e.g.
            0       {GO:0010008, GO:0070062, IPR036865, GO:0048471...
            1       {GO:0006351, GO:0070062, GO:0007623, GO:004851...
            2       {GO:0019888, GO:0006470, GO:0001754, GO:009024...
    label_type : str
        The prefix to be appended (without the "@" separator).

    Returns
    -------
    labeled_annotations : array-like (pandas Series)
        A new pandas Series where all annotations have received a prefix.
    """
    labeled_annotations = merged_annotations_column.map(
        lambda x: set([label_type + '@' + i for i in x]))
    return labeled_annotations


def _create_dummies(interaction_dataframe, columns):
    """Creates a binary array describing the occurrence of pairs. Both axis
    contain annotation labels, 1's signify the existence of a pairing and 0
    the absence.

    NOTE: requires NaNs to be replaced by empty sets!

    Parameters
    ----------
    interaction_dataframe : [type]
        [description]
    columns : [type]
        [description]
    """

    # replace NaNs with empty sets
    for i in columns:
        interaction_dataframe.loc[interaction_dataframe[
            i].isnull(), i] = interaction_dataframe.loc[interaction_dataframe[
                i].isnull(), i].apply(lambda x: set())

    dummies_df = pd.concat(
        [
            interaction_dataframe[i].apply(
                lambda x: pd.Series([1] * len(x), index=x)).fillna(
                    0, downcast='infer') for i in columns
        ],
        axis=1)

    return dummies_df


def find_pairwise(interaction_dataframe,
                  column_A,
                  column_B,
                  propagate=False,
                  go_dict=None):
    """Retrieve pairwise associations between the annotation sets of two
    interacting proteins.

    Each supplied column (one each for the two interacting proteins) should
    contain a set of annotations. These columns should be constructed using the
    merge_annotations() function.

    The hierarchy of the GO terms are taken into account when the propagate
    flag is set to True (default). This means that an association between two
    terms will be considered the same as an association between any child terms
    of the original two terms.

    Parameters
    ----------
    interaction_dataframe : DataFrame
        DataFrame containing protein-protein interactions and annotations.
    columns : list
        A list of annotation columns between which associations should be
        found. Expects the column contents to be sets of annotations:
            0       {p@GO:0019012, p@GO:0019031, p@GO:0055036, p@G...
            1       {p@IPR013672, p@IPR001889, p@IPR027417, p@GO:0...
            2       {p@IPR013672, p@IPR001889, p@IPR027417, p@GO:0...
    propagate : bool, optional
        Whether or not to propagate GO labels.
        (the default is True, which does not propagate the terms)
    go_dict : dict, optional
        A dictionary containing the GO hierarchy. Constructed via the
        obo_tools.importOBO() function in the goscripts package.

    Returns
    -------
    results : dict
        A dictionary containing the pairwise association measures for each
        pair of annotations.

    """

    columns = [column_A, column_B]
    # replace NaNs with empty sets
    for i in columns:
        interaction_dataframe.loc[interaction_dataframe[
            i].isnull(), i] = interaction_dataframe.loc[interaction_dataframe[
                i].isnull(), i].apply(lambda x: set())

    # TODO: optional step: label pathogen vs host terms

    # create all possible pairwise combinations

    # first create all pairwise combinations per ppi
    # this returns a series where each element/row is a list containing
    # pairwise tuples that are sorted, i.e. (GO:0016197, GO:0055036) and
    # (GO:0055036, GO:0016197) willl not co-occurr in a list
    interaction_dataframe['annotation_pairs'] = interaction_dataframe[columns].apply(
        lambda x: set([i for i in itertools.product(x[column_A], x[column_B])]),
        axis=1)
    # results in sets of tuple pairs
    #   0       {(h@GO:0017137, p@GO:0019012), (h@GO:0005515, ...
    #   1       {(h@GO:0008380, p@GO:0071897), (h@GO:0051219, ...
    #   2       {(h@GO:0007525, p@GO:0071897), (h@GO:0061053, ...

    # NOTE: if the annotations are labelled as host/pathogen, the "i"s don't
    #       really need to be sorted explicitly. Instead, order is given by
    #       the order in which the two lists are provided.
    # interaction_dataframe['annotation_pairs'] = interaction_dataframe[columns].apply(lambda x: set([tuple(sorted(i)) for i in itertools.product(x[column_A], x[column_B])]), axis=1)

    # next, these lists are joined and duplicates are removed
    pairs = set().union(*interaction_dataframe['annotation_pairs'])

    # this would create too many non-existing combinations
    # annotation_set = set().union(
    #   *interaction_dataframe[columns].values.ravel('K'))
    # list(itertools.combinations(annotation_set, 2))

    # create a set of annotations for each interaction as a whole (i.e. join
    # labels for the two interacting proteins)
    merged_annotations = merge_annotations(interaction_dataframe, columns)
    # results in sets of annotations, across both interacting proteins
    #       0       {h@GO:0017137, h@GO:2001136, p@GO:0019012, h@G...
    #       1       {h@GO:0038023, h@GO:0005654, h@GO:0000381, p@I...
    #       2       {h@GO:0006470, p@IPR013672, h@GO:0090263, h@GO...

    # create dummies dataframe
    # https://stackoverflow.com/questions/29034928/pandas-convert-a-column-of-list-to-dummies
    # dummies_df = _create_dummies(interaction_dataframe, columns)

    # results_dict = {
    #     pair: {
    #         'jaccard': {},
    #         'pmi': {},
    #         'G': {},
    #         'chi2': {},
    #         'fisher': {},
    #         'phi': {},
    #         'min_count': {}
    #     }
    #     for pair in pairs
    # }

    results_dict = {
        # 'jaccard': {},
        'pmi': {},
        'G': {},
        'chi2': {},
        'fisher': {},
        'phi': {},
        'min_count': {},
        'counts': {},
        'depth': {}
    }

    counter = 0

    import time
    start = time.time()

    if not propagate:

        for pair in pairs:

            counter += 1
            if counter % 10000 == 0:
                print('Processed {} lines out of {}'.format(
                    counter, len(pairs)))

            pair_count, label_one_count_exclusive, label_two_count_exclusive, label_one_count, label_two_count, absent_count, total_count = count_presences(
                pair, interaction_dataframe, 'annotation_pairs',
                merged_annotations)

            # results_dict[pair]['pmi'] = _calc_pmi(pair_count, label_one_count, label_two_count)
            # results_dict[pair]['chi2'] = _calc_chi2(pair_count, label_one_count_exclusive, label_two_count_exclusive, absent_count)
            # results_dict[pair]['G'] = _calc_G(pair_count, label_one_count_exclusive, label_two_count_exclusive, absent_count)
            # results_dict[pair]['fisher'] = _calc_fisher(pair_count, label_one_count_exclusive, label_two_count_exclusive, absent_count)

            results_dict['pmi'][pair] = _calc_pmi(pair_count, label_one_count,
                                                  label_two_count, total_count)
            results_dict['chi2'][pair] = _calc_chi2(
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count)
            results_dict['G'][pair] = _calc_G(
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count)
            results_dict['fisher'][pair] = _calc_fisher(
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count)
            results_dict['min_count'][pair] = min([
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count
            ])

        df_list = [
            pd.DataFrame.from_dict(d, orient='index').add_prefix(name + '_')
            for name, d in results_dict.items()
        ]
        results = pd.concat(df_list, axis=1)

        end = time.time()
        print('Pairwise calculations took {}'.format(end - start))

        return results

    else:
        pairs_dict = propagate_pairs(pairs, go_dict)

        for pair, propagated_pair in pairs_dict.items():

            counter += 1
            if counter % 10000 == 0:
                print('Processed {} lines out of {}'.format(
                    counter, len(pairs)))

            pair_count, label_one_count_exclusive, label_two_count_exclusive, label_one_count, label_two_count, absent_count, total_count = count_presences_propagated(
                propagated_pair, interaction_dataframe, column_A, column_B,
                merged_annotations)

            # print(pair)
            # print(np.array([[pair_count, label_one_count_exclusive],
            #     [label_two_count_exclusive, absent_count]]))
            # print('pair_count, label_one_count_exclusive, label_two_count_exclusive, label_one_count, label_two_count, absent_count, total_count')
            # print(pair_count, label_one_count_exclusive, label_two_count_exclusive, label_one_count, label_two_count, absent_count, total_count)

            results_dict['pmi'][pair] = _calc_pmi(pair_count, label_one_count,
                                                  label_two_count, total_count)
            results_dict['chi2'][pair] = _calc_chi2(
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count)
            results_dict['G'][pair] = _calc_G(
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count)
            results_dict['fisher'][pair] = _calc_fisher(
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count)
            results_dict['min_count'][pair] = min([
                pair_count, label_one_count_exclusive,
                label_two_count_exclusive, absent_count
            ])
            results_dict['counts'][pair] = {
                'pair_count': pair_count,
                'label_one_count_exclusive': label_one_count_exclusive,
                'label_two_count_exclusive': label_two_count_exclusive,
                'label_one_count': label_one_count,
                'label_two_count': label_two_count,
                'absent_count': absent_count,
                'total_count': total_count
            }
            results_dict['depth'][pair] = {
                'protein_A':
                go_dict[pair[0][2:]].depth if 'GO' in pair[0] else None,
                'protein_B':
                go_dict[pair[1][2:]].depth if 'GO' in pair[1] else None
            }

        df_list = [
            pd.DataFrame.from_dict(d, orient='index').add_prefix(name + '_')
            for name, d in results_dict.items()
        ]
        results = pd.concat(df_list, axis=1)

        end = time.time()
        print('Pairwise calculations took {}'.format(end - start))

        return results


def count_presences(pair, interaction_dataframe, pairs_column,
                    merged_annotations):
    """
    something
    """

    #TODO: maps > apply! merge bool arrays > multiple queries!

    # count all ppis where pair occurs
    presence_mask = interaction_dataframe[pairs_column].map(
        lambda x: pair in x)
    pair_count = np.sum(presence_mask)

    # count ppis where only 1 label occurs: P(X | Y') or N10
    label_one_count_exclusive = np.sum(
        merged_annotations.map(lambda x: pair[0] in x and pair[1] not in x))
    label_two_count_exclusive = np.sum(
        merged_annotations.map(lambda x: pair[0] not in x and pair[1] in x))

    # count ppis where 1 label occurs, regardless of other label in pair:
    # P(X) or N1+
    label_one_count = np.sum(merged_annotations.map(lambda x: pair[0] in x))
    label_two_count = np.sum(merged_annotations.map(lambda x: pair[1] in x))

    # count ppis lacking either term: P(X',Y') or N00
    absent_count = np.sum(
        merged_annotations.map(
            lambda x: pair[0] not in x and pair[1] not in x))

    # total count of terms
    # TODO: check?
    total_count = interaction_dataframe.shape[0]

    return pair_count, label_one_count_exclusive, label_two_count_exclusive, label_one_count, label_two_count, absent_count, total_count


def count_presences_propagated(pair_prop, interaction_dataframe, column_A,
                               column_B, merged_annotations):
    """
    columns must be host, pathogen order
    """
    """Comparison of methods

    ALL PAIRS
    %time np.sum([np.sum(interaction_dataframe['annotation_pairs'].apply(lambda x: i in x)) for i in list(itertools.product(pair_prop[0], pair_prop[1]))])
    %time np.sum(interaction_dataframe.apply(lambda x: not x['annotations_A'].isdisjoint(pair_prop[0]) and not x['annotations_B'].isdisjoint(pair_prop[1]), axis=1))
    %time np.sum(merged_annotations.map(lambda x: not x.isdisjoint(pair_prop[0]) and not x.isdisjoint(pair_prop[1])))

    CPU times: user 138 ms, sys: 30 µs, total: 138 ms
    Wall time: 136 ms
    CPU times: user 105 ms, sys: 0 ns, total: 105 ms
    Wall time: 105 ms
    CPU times: user 8.65 ms, sys: 0 ns, total: 8.65 ms
    Wall time: 8.46 ms

    ONLY 1 LABEL

    """

    # count all ppis where any pair occurs
    pair_count = np.sum(
        merged_annotations.map(
            lambda x: not x.isdisjoint(pair_prop[0]) and not x.isdisjoint(pair_prop[1])
        ))

    # count ppis where only 1 label occurs: P(X | Y') or N10
    label_one_count_exclusive = np.sum(interaction_dataframe[column_A].map(
        lambda x: not x.isdisjoint(pair_prop[0])
    ) & interaction_dataframe[column_B].map(
        lambda x: x.isdisjoint(pair_prop[1])))
    label_two_count_exclusive = np.sum(interaction_dataframe[column_A].map(
        lambda x: x.isdisjoint(pair_prop[0])
    ) & interaction_dataframe[column_B].map(
        lambda x: not x.isdisjoint(pair_prop[1])))

    # count ppis where 1 label occurs, regardless of other label in pair:
    # P(X) or N1+
    label_one_count = np.sum(
        merged_annotations.apply(lambda x: not x.isdisjoint(pair_prop[0])))
    label_two_count = np.sum(
        merged_annotations.apply(lambda x: not x.isdisjoint(pair_prop[1])))
    # option to not discern between location of label (i.e. host or pathogen)
    # in case there would be identical labels in host and pathogen proteins
    # label_one_count = np.sum(interaction_dataframe[column_A].map(
    #     lambda x: not x.isdisjoint(pair_prop[0])))
    # label_two_count = np.sum(interaction_dataframe[column_B].map(
    #     lambda x: not x.isdisjoint(pair_prop[1])))

    # count ppis lacking either term: P(X',Y') or N00
    absent_count = np.sum(
        merged_annotations.map(lambda x: x.isdisjoint(pair_prop[0])) &
        merged_annotations.map(lambda x: x.isdisjoint(pair_prop[1])))
    # option to not discern between location of label (i.e. host or pathogen)
    # absent_count = np.sum(interaction_dataframe[column_A].map(
    #     lambda x: x.isdisjoint(pair_prop[0])) & interaction_dataframe[column_B]
    #                       .map(lambda x: x.isdisjoint(pair_prop[1])))

    # total count of terms
    total_count = interaction_dataframe.shape[0]

    return pair_count, label_one_count_exclusive, label_two_count_exclusive, label_one_count, label_two_count, absent_count, total_count


# def _jaccard(a, b, c):
#     return a / (a + b + c)


def _calc_pmi(n11, n1plus, nplus1, nplusplus):
    return math.log(
        (n11 / nplusplus) / ((n1plus / nplusplus) * (nplus1 / nplusplus)), 2)


def _calc_chi2(n11, n10, n01, n00):
    # print(n11, n10, n01, n00)
    # print(type(n11))
    contingency_table = np.array([[n11, n01], [n10, n00]])

    # if 0 in contingency_table:
    #     print(contingency_table)
    #     # import ipdb
    #     # ipdb.set_trace()

    chi2, p, df, ex = scipy.stats.chi2_contingency(
        contingency_table, correction=True, lambda_="pearson")
    d = {'chi2': chi2, 'p-value': p, 'df': df, 'exp': ex}
    return d


def _calc_G(n11, n10, n01, n00):
    contingency_table = np.array([[n11, n01], [n10, n00]])
    G, p, df, ex = scipy.stats.chi2_contingency(
        contingency_table, correction=False, lambda_="log-likelihood")
    d = {'G': G, 'p-value': p, 'df': df, 'exp': ex}
    return d


def _calc_fisher(n11, n10, n01, n00):
    contingency_table = np.array([[n11, n01], [n10, n00]])
    oddsratio, p = scipy.stats.fisher_exact(contingency_table)
    d = {'oddsration': oddsratio, 'p-value': p}
    return d


def multiple_testing_correction(results_dataframe,
                                columns,
                                threshold=0.05,
                                method='fdr_bh'):
    """Performs multiple testing correction on a list of p-values.

    Parameters
    ----------
    results_dataframe : DataFrame
        A pandas DataFrame containing one or more columns with p-values.
    columns : list
        A list containing the names of the columns with p-values.

    Returns
    -------
    None
        Modifies the DataFrame in-place by adding a column with corrected
        p-values. The name is constructed by taking the original column name
        and appending the chosen method, e.g. "fisher_p-value_fdr_bh".
    """

    for i in columns:
        p_values = results_dataframe[i].values
        try:
            corrected_p = statsmodels.stats.multitest.multipletests(
                p_values, alpha=threshold, method=method)
        except ValueError:
            print(
                'ERROR: Invalid method for multiple testing correction. Accepted options include: fdr_bh, bonferroni and any others defined by statsmodels.stats.multitest.multipletests().'
            )
            sys.exit(1)

        results_dataframe[i + '_' + method] = corrected_p[1]


def propagate_pairs(pairs, go_dict):
    """Propagates pairs of annotations.

    For each pair in an array of annotation pairs, the GO annotations will
    be replaced by a set of their (recursive) child terms (including itself).

    Other types of annotations are left untouched, but converted to a 1-member
    set.

    Parameters
    ----------
    pairs : array
        An array of sorted tuples of annotations, one for the host and one for
        the pathogen, e.g. ('h@GO:0030133', 'p@IPR009304')
    go_dict : dict
        A dictionary containing the GO hierarchy. Constructed via the
        obo_tools.importOBO() function in the goscripts package.

    Returns
    -------
    propagated_pairs : dict
        A dictionary of sorted tuples of annotation sets, one for the host and
        one for the pathogen. Each element in the tuple consists of a set of
        terms, e.g. a GO term and all of its descendants.
    """

    propagated_pairs = {}
    for pair in list(pairs):
        new_pair = []
        for term in pair:
            if 'GO' in term:
                prefix = term[:2]
                go_object = go_dict.get(term[2:])
                if None:
                    new_pair.append(term)
                else:
                    child_terms = [
                        prefix + i for i in go_object.recursive_children
                    ]
                    propagation_set = set(child_terms) | set([term])
                    new_pair.append(propagation_set)
            else:
                new_pair.append({term})
        propagated_pairs[pair] = tuple(new_pair)

    return propagated_pairs


"""Parallellization attempts

short_pairs = list(pairs)[0:1000]

list_pairs = list(pairs)

from math import sqrt

from joblib import Parallel, delayed

%time Parallel(n_jobs=7)(delayed(count_stuff)(i, ppi_df, columns, merged_annotations) for i in short_pairs)

%time Parallel(n_jobs=7,backend="threading")(delayed(count_stuff)(i, ppi_df, columns, merged_annotations) for i in short_pairs)

%time [count_stuff(i, ppi_df, columns, merged_annotations) for i in short_pairs]

CPU times: user 13min 15s, sys: 1min 3s, total: 14min 19s
Wall time: 18min 4s
CPU times: user 17.9 s, sys: 441 ms, total: 18.4 s
Wall time: 17 s
CPU times: user 13.6 s, sys: 5 ms, total: 13.6 s
Wall time: 13.6 s

"""
