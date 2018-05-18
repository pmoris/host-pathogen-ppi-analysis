#!/usr/bin/env python3
"""
Module to parse NCBI taxon ID database.

Can also be run as a stand-alone script to generate all child taxon IDs by
providing the taxdump directory and desired parent taxid as arguments.

Usage: python3 src/dataprep/taxonid.py data/raw/taxdump/ 10292 output-file
"""

import os
import sys
from pathlib import Path


def parse_taxid_names(file_path):
    """
    Parse the names.dmp file and output a dictionary mapping names to taxids
    (multiple different keys) and taxids to scientific names.

    Parameters
    ----------
    file_path : str
        The path to the names.dmp file.

    Returns
    -------
    name2taxid : dict
        Keys are all possible names and values are taxids.
    taxid2name : dict
        Keys are taxids and values are scientific names.
    """

    names = Path(file_path)

    with names.open() as f:
        lines_processed = 0
        name2taxid = {}
        taxid2name = {}
        for line in f:
            lines_processed += 1
            if lines_processed % 1000000 == 0:
                print('processing line', str(lines_processed))
            entries = [entry.strip() for entry in line.split('|')]
            name2taxid[entries[1]] = entries[0]
            if 'scientific name' in line:
                taxid2name[entries[0]] = entries[1]

    return name2taxid, taxid2name


def parse_taxid_nodes(file_path):
    nodes = Path(file_path)

    with nodes.open() as f:
        lines_processed = 0
        taxid2parent = {}
        taxid2rank = {}
        for line in f:
            lines_processed += 1
            if lines_processed % 1000000 == 0:
                print('processing line', str(lines_processed))
            entries = [entry.strip() for entry in line.split('|')]
            taxid2parent[entries[0]] = entries[1]
            taxid2rank[entries[0]] = entries[2]

    return taxid2parent, taxid2rank


def name_search(partial_name, name2taxid):
    return [
        taxid for key, taxid in name2taxid.items()
        if partial_name.lower() in key.lower()
    ]


def find_rank(species, rank):
    return find_rank_recursive(name2taxid[species], rank)


def find_rank_recursive(taxid, rank):
    # if the taxid in the argument is the order, then return the name
    if taxid2rank[taxid] == rank:
        return taxid2name[taxid]
    # otherwise, run the same function with the parent of the taxid
    else:
        return find_rank_recursive(taxid2parent[taxid], rank)


def retrieve_parents(taxid):
    return retrieve_parents_recursive(taxid, [])


def retrieve_parents_recursive(taxid, parent_list):
    if taxid == 1:
        return parent_list
    else:
        parent = taxid2parent[taxid]
        parent_list.append(parent)
        return retrieve_parents_recursive(parent, parent_list)


def find_lca(taxid1, taxid2):
    parents1 = retrieve_parents(taxid1)
    parents2 = retrieve_parents(taxid2)
    for parent in parents1:
        if parent in parents2:
            return parent


def find_lca_group(taxid_list):
    parent_list = [retrieve_parents(taxid) for taxid in taxid_list]
    # first loop ensures order...
    for i in parent_list[0]:
        if all(i in j for j in parent_list[1:]):
            return (i)


def create_parent2child_dict(taxid2parent_dict):
    parent2child = {}
    for taxid, parent in taxid2parent_dict.items():
        if parent not in parent2child.keys():
            parent2child[parent] = []
        parent2child[parent].append(taxid)
    return parent2child


def get_children(taxid, parent2child):
    # get the children of the current node or empty list if it doesn't exist
    # requires copy of list to stop appending of childs in original dictionary
    children = parent2child.get(taxid, [])[:]
    # add all the children of this node to the list
    # then add all of their children as well
    # requires copy of list otherwise it will grow and the iteration will create duplicates
    for child in children[:]:
        children.extend(get_children(child, parent2child))
    return children


# def get_children(taxid, parent2child):
#     print('Retrieving taxids falling under',taxid,taxid2name[taxid])
#     # get the children of the current node or empty list if it doesn't exist
#     children = parent2child.get(taxid, [])
#     # initiate set to store taxids and
#     # add the children to it
#     result_set = set(children)
#     # Recursively loop through children of initial taxid and add them to result_set.
#     for child in children:
#         get_children_recursive(child, parent2child, result_set)
#     return result_set
#
# def get_children_recursive(taxid, parent2child, result_set):
#     # get the children of the current node or empty list if it doesn't exist
#     children = parent2child.get(taxid, [])
#     # Add them to result_set
#     result_set.update(children)
#     # Recursively call function
#     for child in children:
#         get_children_recursive(child, parent2child, result_set)
#


def write_taxids(taxid_list, taxid2name_dict, out_file):
    with open(Path(out_file), 'w') as out:
        for i in taxid_list:
            out.write(str(i) + '|' + taxid2name_dict[i] + '\n')


# if __name__ == "__main__":
#     '''
#     Generate all child taxon IDs by providing the taxdump directory and
#     desired parent taxid as arguments.
#     '''
#     try:
#         taxdump_dir = Path(sys.argv[1])
#     except IndexError:
#         print('Incorrect path provided.')
#     else:
#         print('Parsing taxdump files...')
#         names = taxdump_dir / 'names.dmp'
#         name2taxid, taxid2name = parse_taxid_names(str(names))
#         nodes = taxdump_dir / 'nodes.dmp'
#         taxid2parent, taxid2rank = parse_taxid_nodes(str(nodes))
#         parent2child = create_parent2child_dict(taxid2parent)

#         taxid = sys.argv[2]
#         print('Retrieving all child taxa of', taxid)
#         children = get_children(taxid, parent2child)
#         if len(sys.argv) > 3:
#             out_path = sys.argv[3]
#         else:
#             out_name = 'child_taxids_of_' + str(taxid) + '.txt'
#             out_path = Path(os.path.abspath(__file__)).parents[
#                 2] / 'data/interim' / out_name
#             print('Saving output as', out_path)
#         save_children(children, str(out_path)) #changed to write_taxids()
#         print('Saved output to', str(out_path))

        # except IndexError:
        #     print('No taxid was provided.')
        # except ValueError:
        #     print('Taxids should be provided as integers')
        # except LookupError:
        #     print('Taxid was not found.')

        # #
        # parent2child = {taxid: get_children(taxid) for taxid in
        #                 [6278, 6296, 6295, 6274, 119089, 6231, 1206794, 33317, 33213, 6072, 33208, 33154, 2759, 131567,
        #                  1]}
        #
        # # print(parent2child[6274])

    #
    #
    # print(find_rank('Panagrolobus vanmegenae', 'family'))
    #
    #
    #
    # print(retrieve_parents(6279))
    #
    #
    #
    # print('lca6279,6239',find_lca(6279,6239))
    # # print('lca,6239,6279',find_lca(6239,6279))
    # # print('find_lca(6279,find_lca(6239,135651))',find_lca(6279,find_lca(6239,135651)))
    # # print('find_lca(6239,find_lca(6279,135651))',find_lca(6239,find_lca(6279,135651)))
    # # print('find_lca(135651,find_lca(6279,6239))',find_lca(135651,find_lca(6279,6239)))
    # # print('find_lca(135651,find_lca(6239,6279))',find_lca(135651,find_lca(6239,6279)))
    #
    #
    #
    # print(find_lca_group([6239,6279,135651]))
    # print(retrieve_parents(10090))
    # print(find_lca_group([6239,6279,135651]) in retrieve_parents(10090))

    #
    # print('testing')
    # print(taxid2parent)
    # if 10296 in taxid2parent.items():
    #     print('yes')
    #
    # parent2child = create_parent2child_dict(taxid2parent)
    #
    #     # print(parent2child[parent])
    #
    # print(taxid2parent[1982744])
    # print(taxid2parent[1982743])
    # print('parent2child', parent2child[35247])
    #
    #
    # print('p2c',parent2child[10292])
    #
    # print('taxid2p',taxid2parent[10296])
    # # print('pare2chil',parent2child[10296])
    #
    # print('get',get_children(10292))
    #
    # # parent2child = {taxid: get_children(taxid) for taxid in [6278, 6296, 6295, 6274, 119089, 6231, 1206794, 33317, 33213, 6072, 33208, 33154, 2759, 131567, 1]}
    #
    # # print(parent2child[6274])
