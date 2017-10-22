#!/usr/bin/env python3

"""
Script to filter relevant data out of the complete InterPro data set.
"""

import os
import sys
from pathlib import Path

sys.path.append('..')

import pandas as pd

import ppi_tools

# https://stackoverflow.com/questions/8951255/import-script-from-a-parent-directory
# http://python-notes.curiousefficiency.org/en/latest/python_concepts/import_traps.html

repo_source_path = Path(os.path.abspath(__file__)).parents[2]

df_virhost = ppi_tools.ppi_import.read_mitab_virhost(repo_source_path / 'data/raw/ppi_data/VirHostNet_January2017.txt')
    # r'../ppi_data/VirHostNet_January_2017.txt')

df_hpidb2 = ppi_tools.ppi_import.read_mitab_hpidb2(repo_source_path / 'data/raw/ppi_data/ppi_data/hpidb2_March14_2017_mitab_plus.txt')
    # r'../ppi_data/hpidb2_March14_2017_mitab_plus.txt')

df_phisto = ppi_tools.ppi_import.read_mitab_phisto(repo_source_path / 'data/raw/ppi_data/phisto_Jan19_2017.csv',
                                        repo_source_path / 'data/raw/ppi_data/mi.obo')
    # r'../ppi_data/phisto_Jan19_2017.csv',
    #                           r'../ppi_data/mi.obo')

df_concat = ppi_tools.main.concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto])

# unique_ac = pd.unique(df_concat[['xref_A', 'xref_B']].values.ravel())

unique_ac = set(pd.unique(df_concat['xref_B'].str.extract('^.*:(\w*)-?',
                                             expand=False).append(df_concat['xref_A'].str.extract('^.*:(\w*)-?',
                                                                                                  expand=False))))

protein2ipr_file = repo_source_path / 'data/raw/interpro_data/protein2ipr.txt'
    # r'../domain-motif/protein2ipr.txt')

protein2ipr_reduced_file = repo_source_path / 'data/interim/interpro_data/protein2ipr_filtered2.txt'
    # r'../domain-motif/protein2ipr_filtered.txt')

uniprot2interpro_dict = {}

with protein2ipr_file.open() as protein2ipr:
    with protein2ipr_reduced_file.open('w') as protein2ipr_reduced:
        outer_counter = 0
        inner_counter = 0
        for line in protein2ipr:
            if outer_counter % 10000000 == 0:
                print('Processed',outer_counter,'lines')
            if line.split('\t')[0] in unique_ac:
                inner_counter += 1
                protein2ipr_reduced.write(line)
                if inner_counter % 1000 == 0:
                    print('Written',inner_counter,'lines.')
            outer_counter += 1
