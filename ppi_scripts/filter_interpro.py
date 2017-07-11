#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

import parse_ppi

df_virhost = parse_ppi.read_mitab_virhost(r'../ppi_data/VirHostNet_January_2017.txt')

df_hpidb2 = parse_ppi.read_mitab_hpidb2(r'../ppi_data/hpidb2_March14_2017_mitab_plus.txt')

df_phisto = parse_ppi.read_mitab_phisto(r'../ppi_data/phisto_Jan19_2017.csv',
                              r'../ppi_data/mi.obo')

df_concat = parse_ppi.concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto])

# unique_ac = pd.unique(df_concat[['xref_A', 'xref_B']].values.ravel())

unique_ac = set(pd.unique(df_concat['xref_B'].str.extract('^.*:(\w*)-?',
                                             expand=False).append(df_concat['xref_A'].str.extract('^.*:(\w*)-?',
                                                                                                  expand=False))))

protein2ipr_file = Path(r'../domain-motif/protein2ipr.txt')

protein2ipr_reduced_file = Path(r'../domain-motif/protein2ipr_filtered.txt')

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
