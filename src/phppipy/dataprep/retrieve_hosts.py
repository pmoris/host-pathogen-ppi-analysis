# import pandas as pd

# from ppi_tools import id_mapper
# import ..ppi_tools.label_go
# import ppi_tools.label_interpro
# import ppi_tools.import

# from data_prep import retrieve_taxids

# import goscripts

# from ppi_tools import import

# works if called from script directory
# import sys
# sys.path.append('..')

# works with python -m dataprep.quick-filter call from src directory
from ppi_tools import import

# Read in PPI datasets

# glob folder

df_virhost = import.read_mitab_virhost('../data/raw/ppi_data/VirHostNet_January_2017.txt')
# df_hpidb2 = import.read_psi_mi_tab('data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt', 'hpidb2')
# df_phisto = import.read_mitab_phisto('data/raw/ppi_data/phisto_Jan19_2017.csv',
#                                             'data/raw/ppi_data/mi.obo')
# df_intact = import.read_psi_mi_tab('data/raw/ppi_data/intact_virus_2017_12_1.txt', 'intact')

# # Concatenate the different sources
# print('Concatenating PPI datasets...')
# df_concat = concat_interaction_datasets([df_hpidb2, df_virhost, df_phisto, df_intact])

# # Remove duplicate interaction pairs (including different detection methods and publications)
# # https://stackoverflow.com/a/41650846
# # https://stackoverflow.com/questions/33042777/
# # Note that this will result in the first dataset (e.g. hpidb2) having priority over the others.
# df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted'], keep='first')
# # df_concat_dedup = df_concat.drop_duplicates(subset=['xref_partners_sorted', 'taxid_B', 'taxid_A'], keep='first')
# df_concat_dedup = df_concat_dedup.reset_index(drop=True)


# # Retrieve Herpesviridae (taxid:10292) list, see retrieve_taxids.py script to generate child taxids
# # TODO: import from retrieve_Taxids and create on the spot
# # TODO: then combine this code into a function to subset a given dataframe for a given taxid and its children
# # TODO: and a host list or default to all hosts
# # TODO: note that this filtering will generally result in intra-host interactions being omitted,
# # TODO: while retaining intra-viral ones
# with Path('data/interim/child_taxids_of_10292.txt').open() as taxid_file:
#     herpes_taxids = [str('taxid:' + line.split('|')[0]) for line in taxid_file]
# print(r'Retrieving Herpes taxids from data/interim/child_taxids_of_10292.txt')