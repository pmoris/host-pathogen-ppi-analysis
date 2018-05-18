from ..phppipy.ppi_tools import import

# with open('data/raw/ppi_data/hpidb2_March14_2017_mitab.txt', encoding='ISO-8859-1') as f:
#     first_line = f.readline()
#     print(first_line)

# df = import.read_mi_tab('data/raw/ppi_data/VirHostNet_January_2017.txt')
# print(df.head())
# df = import.read_mi_tab('data/raw/ppi_data/hpidb2_March14_2017_mitab.txt')
# print(df.head())
df = import.read_mi_tab('data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt')
print(df.head())

print(df.columns.values)
import sys
print(__name__)
print(__package__)
import os
# print(os.path)