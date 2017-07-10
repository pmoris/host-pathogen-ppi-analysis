def read_mitab_hpidb2_manual(filepath):
    file = Path(filepath)
    with file.open(encoding="ISO-8859-1") as f:
        dict_list = []
        for line in f:
            l = line.strip().split('\t')
            keys = ['protein_xref_1', "protein_xref_2", "alternative_identifiers_1", "alternative_identifiers_2",
                    "protein_alias_1", "protein_alias_2", "detection_method", "author_name", "pmid", "protein_taxid_1",
                    "protein_taxid_2", "interaction_type", "source_database_id", "database_identifier", "confidence",
                    "protein_xref_1_unique", "protein_xref_2_unique", "protein_taxid_1_cat", "protein_taxid_2_cat",
                    "protein_taxid_1_name", "protein_taxid_2_name", "protein_seq1", "protein_seq2", "source_database",
                    "protein_xref_1_display_id", "protein_xerf_2_display_id"]
            dict_list.append({k: v for k, v in zip(keys, l)})
    df = pd.DataFrame(dict_list)
    return df


def read_mitab_pandas_chunk(filepath):
    file = Path(filepath)
    df = pd.read_csv(file, sep='\t', encoding="ISO-8859-1", iterator=True, chunksize=1000)
    df_concat = pd.concat([chunk for chunk in df])
    return df_concat

'''
# Dask implementation
import dask.dataframe as dd
import dask.multiprocessing
def read_mitab_dask(filepath):
    file = Path(filepath)
    df = dd.read_csv(file, sep='\t', blocksize=1000000)
    df = df.compute(get=dask.multiprocessing.get)
    return df

s = read_mitab(r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi-data/'
               r'hpidb2_March14_2017_mitab_plus.txt')
'''

# Add unique identifier for interaction pairs
xref_partners_sorted_array = np.sort(np.stack((df_concat.xref_A, df_concat.xref_B), axis=1), axis=1)
df_concat['xref_partners_sorted'] = list(map(tuple, xref_partners_sorted_array))
df_concat['xref_partners_sorted'] = pd.Series(map(tuple, xref_partners_sorted_array))

# sorting the following does not work:
# df_concat['xref_partners_sorted'] = list(zip(df_concat.xref_A, df_concat.xref_B))

# slower alternative, returns series of lists
# f_concat.apply(lambda x: sorted(x[['xref_A', 'xref_B']]), axis=1)
# other options
# https://stackoverflow.com/questions/40187874/python-pandas-two-columns-with-same-values-alphabetically-sorted-and-stored
# use df[['col1','col2']].values to get array and use np.sort
# or add .tolist() to get list and use sorted(i) for i in list
# or sort in-place df.values.sort(axis=1)
# pandas.DataFrame.sort_values




# reorder pathogen host columns

# Slower alternative
# lambdafunc = lambda x: pd.Series([x['xref_B'],x['xref_A']])
# df_herpes.loc[host_position_mask, ['xref_A', 'xref_B']] = df_herpes.loc[host_position_mask].apply(lambdafunc, axis=1)




# reordering printouts
'''
# print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_A.values)
# print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_B.values)
#
# print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_A').size())
# print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_B').size())
#
# print(df_herpes.groupby('taxid_A').size())
# print(df_herpes.groupby('taxid_B').size())
'''

print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_A.values)
print('taxid:9606' in df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].taxid_B.values)

print('\n\n\n\n\nvirhosttaxidssizesgrouped\n\n\n\n')
print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_A').size())
print(df_herpes.loc[df_herpes['origin'] == 'VirHostNet2'].groupby('taxid_B').size())






# Equivalent to label_host_pathogen_origin() functoin that loops through
# column-taxid zip and redefines a new lambda function each time.

 def label_row_via_apply(row, pathogen_set, GO_columns, taxid_columns):
    labelled_columns = pd.Series()
    for i, j in zip(GO_columns, taxid_columns):
        labelled_columns = labelled_columns.append(pd.Series(label_term(row[i], row[j], pathogen_set)))
    return labelled_columns

def label_host_pathogen_origin(interaction_dataframe, pathogen_set, GO_columns=list(('xref_A_GO', 'xref_B_GO')),
                               taxid_columns=list(('taxid_A', 'taxid_B'))):
    interaction_dataframe[GO_columns] = interaction_dataframe.apply(label_row_via_apply,
                                                                    args=(pathogen_set, GO_columns, taxid_columns),
                                                                    axis=1)

    # original code
    '''
    lambda_go_label = lambda x: pd.Series([label_GO_term(x['xref_A_GO'], x['taxid_A'], host_taxids),
                                           label_GO_term(x['xref_B_GO'], x['taxid_B'], host_taxids)])

    df_herpes[['xref_A_GO', 'xref_B_GO']] = df_herpes.apply(lambda_go_label, axis=1)
    '''



    # # Check number of within and between interactions
    # def check_intra_species_interactions(interaction_dataframe):
    #     print('Checking for intra-species interactions')
    #     print(np.where(df_hpidb2.taxid_A == df_hpidb2.taxid_B))  # none
    #     print(np.where(df_virhost.taxid_A == df_virhost.taxid_B))  # many
    #     print(np.where(df_virhost.taxid_A == df_virhost.taxid_B)[0].size)  # 28664 / 752190
    #     print(df_virhost.loc[(df_virhost.taxid_A == df_virhost.taxid_B), ['taxid_A', 'taxid_B']].shape)
    #






# pathogen to higher group mapping


    # # ateline_gh3 = ['85618']
    # pathogen_group_dict['bovine_ah1'] = ['10320'] + retrieve_taxids.get_children('10320',parent2child)
    # # bovine_gh4 = ['10385'].extend(retrieve_taxids.get_children('10385',parent2child))
    # pathogen_group_dict['bovine_hv1'] = ['79889'] + retrieve_taxids.get_children('79889', parent2child)
    # # Eleph_hv1 = ['654902'] # actual parent taxa = 146015
    # pathogen_group_dict['epstein_barr'] = ['10376'] + retrieve_taxids.get_children('10376',parent2child)  # or human gammaherpesvirus 4
    # pathogen_group_dict['equid_av1'] = ['10326'] + retrieve_taxids.get_children('10326', parent2child)  # or equine hv1
    # pathogen_group_dict['equid_gv2'] = ['12657'] + retrieve_taxids.get_children('12657', parent2child)  # or equine hv2
    # pathogen_group_dict['gallid_av2'] = ['10390'] + retrieve_taxids.get_children('10390', parent2child)
    # pathogen_group_dict['human_hsv1'] = ['10298'] + retrieve_taxids.get_children('10298', parent2child)  # or human av1
    # pathogen_group_dict['saimiri_gv2'] = ['10381'] + retrieve_taxids.get_children('10381',parent2child)  # or herpesvirus saimiri
    # pathogen_group_dict['human_av2'] = ['10310'] + retrieve_taxids.get_children('10310', parent2child)
    # pathogen_group_dict['human_av3'] = ['10335'] + retrieve_taxids.get_children('10335',parent2child)  # or varicella-zoster
    # pathogen_group_dict['human_bv5'] = ['10359'] + retrieve_taxids.get_children('10359',parent2child)  # or human cytomegalovirus
    # pathogen_group_dict['human_gv8'] = ['37296'] + retrieve_taxids.get_children('37296',parent2child)  # also contains ape viruses?
    # pathogen_group_dict['human_bv6A'] = ['32603'] + retrieve_taxids.get_children('32603', parent2child)
    # pathogen_group_dict['human_bv6B'] = ['32604'] + retrieve_taxids.get_children('32604', parent2child)
    # pathogen_group_dict['murid_bv1'] = ['10366'] + retrieve_taxids.get_children('10366', parent2child)
    # pathogen_group_dict['murid_gv4'] = ['33708'] + retrieve_taxids.get_children('33708',parent2child)  # or murine hv 68
    # pathogen_group_dict['papiine_gv1'] = ['106332'] + retrieve_taxids.get_children('106332', parent2child)
    # pathogen_group_dict['suid_av1'] = ['10345'] + retrieve_taxids.get_children('10345',parent2child)  # or pseudorabies virus


'''
    print(sorted([taxid2name[i.split(':')[1]] for i in np.sort(np.setdiff1d(all_taxids, host_taxids))]))


pd.Series(np.sort(df_herpes['pathogen_type'].unique())).apply(lambda x: [x, name2taxid[x]])    
\
0                                  [Ateline gammaherpesvirus 3, 85618]
1                                   [Bovine alphaherpesvirus 1, 10320]
2                                   [Bovine gammaherpesvirus 4, 10385]
3                                 [Bovine herpesvirus type 1.1, 79889]
4                 [Bovine herpesvirus type 1.1 (strain Cooper), 10323]
5                   [Bovine herpesvirus type 1.1 (strain Jura), 31518]
6                   [Bovine herpesvirus type 1.1 (strain P8-2), 10324]
7                    [Elephantid herpesvirus 1 (isolate Kiba), 654902]
8                             [Epstein-barr virus strain ag876, 82830]
9                                    [Equid alphaherpesvirus 1, 10326]
10                      [Equid herpesvirus type 2 strain 86/87, 82831]
11                    [Equine herpesvirus type 1 (strain AB4P), 31520]
12              [Equine herpesvirus type 1 (strain Kentucky A), 10329]
13                                  [Gallid alphaherpesvirus 2, 10390]
14                  [Herpes simplex virus (type 1 / strain 17), 10299]
15                   [Herpes simplex virus (type 1 / strain F), 10304]
16              [Herpes simplex virus (type 1 / strain Patton), 10308]
17                            [Herpesvirus saimiri (strain 11), 10383]
18                           [Herpesvirus saimiri (strain 488), 10384]
19                                   [Human alphaherpesvirus 1, 10298]
20                                   [Human alphaherpesvirus 2, 10310]
21                                   [Human alphaherpesvirus 3, 10335]
22                                    [Human betaherpesvirus 5, 10359]
23                                   [Human gammaherpesvirus 4, 10376]
24                                   [Human gammaherpesvirus 8, 37296]
25                             [Human herpesvirus 1 strain KOS, 10306]
26                             [Human herpesvirus 2 strain 333, 10313]
27                            [Human herpesvirus 2 strain HG52, 10315]
28                           [Human herpesvirus 3 strain Dumas, 10338]
29                    [Human herpesvirus 3 strain Oka vaccine, 341980]
30                           [Human herpesvirus 4 strain B95-8, 10377]
31                           [Human herpesvirus 5 strain AD169, 10360]
32                         [Human herpesvirus 5 strain Merlin, 295027]
33                           [Human herpesvirus 5 strain Towne, 10363]

34                            [Human herpesvirus 6 (strain GS), 10369]
35                   [Human herpesvirus 6 (strain Uganda-1102), 10370]
36                             [Human herpesvirus 6 strain Z29, 36351]

37                           [Human herpesvirus 8 strain GK18, 868565]
38                                [Human herpesvirus 8 type M, 435895]

39              [Marek's disease herpesvirus type 1 strain MD5, 10389]

40                                    [Murid betaherpesvirus 1, 10366]

41                                   [Murid gammaherpesvirus 4, 33708]

42                       [Murine cytomegalovirus (strain K181), 69156]
43                      [Murine cytomegalovirus (strain Smith), 10367]

44                                [Papiine gammaherpesvirus 1, 106332]

45                                    [Suid alphaherpesvirus 1, 10345]
46    [Suid herpesvirus 1 (strain Indiana-Funkhauser / Becker), 31523]
dtype: object
(8124, 35)


    print(df_herpes.groupby('pathogen_type').size())

    pathogen_type
Ateline gammaherpesvirus 3                                    1
Bovine alphaherpesvirus 1                                    13
Bovine gammaherpesvirus 4                                     1
Bovine herpesvirus type 1.1                                   6
Bovine herpesvirus type 1.1 (strain Cooper)                   6
Bovine herpesvirus type 1.1 (strain Jura)                     5
Bovine herpesvirus type 1.1 (strain P8-2)                     3
Elephantid herpesvirus 1 (isolate Kiba)                       1
Epstein-barr virus strain ag876                            2160
Equid alphaherpesvirus 1                                      8
Equid herpesvirus type 2 strain 86/87                         4
Equine herpesvirus type 1 (strain AB4P)                       2
Equine herpesvirus type 1 (strain Kentucky A)                 6
Gallid alphaherpesvirus 2                                     4
Herpes simplex virus (type 1 / strain 17)                   925
Herpes simplex virus (type 1 / strain F)                      4
Herpes simplex virus (type 1 / strain Patton)                 1
Herpesvirus saimiri (strain 11)                               5
Herpesvirus saimiri (strain 488)                              4
Human alphaherpesvirus 1                                     34
Human alphaherpesvirus 2                                      3
Human alphaherpesvirus 3                                    133
Human betaherpesvirus 5                                       6
Human gammaherpesvirus 4                                    318
Human gammaherpesvirus 8                                    235
Human herpesvirus 1 strain KOS                                1
Human herpesvirus 2 strain 333                                8
Human herpesvirus 2 strain HG52                              14
Human herpesvirus 3 strain Dumas                              4
Human herpesvirus 3 strain Oka vaccine                       48
Human herpesvirus 4 strain B95-8                           2420
Human herpesvirus 5 strain AD169                             64
Human herpesvirus 5 strain Merlin                           117
Human herpesvirus 5 strain Towne                             21
Human herpesvirus 6 (strain GS)                               1
Human herpesvirus 6 (strain Uganda-1102)                      2
Human herpesvirus 6 strain Z29                                3
Human herpesvirus 8 strain GK18                             369
Human herpesvirus 8 type M                                  281
Marek's disease herpesvirus type 1 strain MD5                17
Murid betaherpesvirus 1                                     330
Murid gammaherpesvirus 4                                    432
Murine cytomegalovirus (strain K181)                         34
Murine cytomegalovirus (strain Smith)                         3
Papiine gammaherpesvirus 1                                    4
Suid alphaherpesvirus 1                                      59
Suid herpesvirus 1 (strain Indiana-Funkhauser / Becker)       4
dtype: int64
(8124, 35)

    '''





def check_host_location(interaction_dataframe):
    print('Checking for presence of humans in taxid A and B columns')
    print('taxid:9606' in df_hpidb2.taxid_A.values)  # true
    # values returns array, faster
    # unique is array, but needs to uniqueify firt
    # print('taxid:9606' in df_hpidb2.taxid_A.unique())
    # print(df_hppid2.taxid_A.isin(['taxid:9606']))
    # print(df_hpidb2.taxid_A.str.contains('9606').any())
    # print(df_hpidb2.taxid_A.str.contains('9606').sum())
    # print(np.any(['9606' in taxid for taxid in df_hpidb2.taxid_A.unique()]))
    # print(np.any(['9606' in taxid for taxid in df_hpidb2.taxid_A.values]))
    # print(len(df_hpidb2.taxid_A))
    print('taxid:9606' in df_hpidb2.taxid_B.unique())  # false
    print('taxid:9606' in df_virhost.taxid_A.unique())  # true
    print('taxid:9606' in df_virhost.taxid_B.unique())  # true


def count_taxids(interaction_dataframe):
    print('Checking number of taxids in A and B')
    print(df_hpidb2.groupby('taxid_A').size())  # contains non-human hosts
    print(df_hpidb2.groupby('taxid_B').size())
    print(df_virhost.groupby('taxid_A').size())
    print(df_virhost.groupby('taxid_B').size())


def check_intra_human_interactions(interaction_dataframe):
    print('Checking for human-human interactions')
    print(df_hpidb2.loc[(df_hpidb2.taxid_A.str.contains('9606')) & (df_hpidb2.taxid_B.str.contains('9606')),
                        ['taxid_A', 'taxid_B']].shape)  # none
    print(df_virhost.loc[(df_virhost.taxid_A.str.contains('9606')) & (df_virhost.taxid_B.str.contains('9606')),
                         ['taxid_A', 'taxid_B']].shape)  # 26732 / 28664


print('HPIDB2 A proteins in VirHost A list')
a = df_hpidb2.loc[df_hpidb2.xref_A.isin(df_virhost.xref_A), ['xref_A']].xref_A.unique()
print('HPIDB2 A proteins in VirHost B list')
b = df_hpidb2.loc[df_hpidb2.xref_A.isin(df_virhost.xref_B), 'xref_A'].unique()
print('HPIDB2 A proteins in VirHost A or B list')
print(np.unique(np.append(a, b)).shape)
print('HPIDB2 A proteins in VirHost A or B list')
print(df_hpidb2.loc[(df_hpidb2.xref_A.isin(df_virhost.xref_A) | df_hpidb2.xref_A.isin(df_virhost.xref_B)),
      :].xref_A.unique().shape)
print('HPIDB2 A proteins unique count')
print(df_hpidb2.xref_A.unique().shape)
print('HPIDB2 A proteins not in VirHost A and not in B list (not in (A or B))')
print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A) | df_hpidb2.xref_A.isin(df_virhost.xref_B)),
      :].xref_A.unique().shape)
print('HPIDB2 B proteins not in VirHost A and not in B list')
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A) | df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins not in VirHost A or not in B list (')
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) | ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins not in VirHost A and not in B list (not in (A or B))')
print(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins in VirHost A or in B list')
print(df_hpidb2.loc[(df_hpidb2.xref_B.isin(df_virhost.xref_A)) | (df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 B proteins in VirHost A or in B list')
print(df_hpidb2.loc[(df_hpidb2.xref_B.isin(df_virhost.xref_A) | df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print('HPIDB2 A proteins not in VirHost A and not in B list OR B proteins not in VirHost A and not in B list')
# entries in hpidb2 in either A or B that don't have any corresponding entry in virhost
print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_B.unique().shape)
print(df_hpidb2.loc[(~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B))) |
                    (~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B))),
      :].xref_B.unique().shape)

print(df_hpidb2.loc[~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)) |
                    ~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)),
      :].xref_A.unique().shape)

print(
    np.unique(df_hpidb2.loc[~(df_hpidb2.xref_B.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_B.isin(df_virhost.xref_B)) |
                            ~(df_hpidb2.xref_A.isin(df_virhost.xref_A)) & ~(df_hpidb2.xref_A.isin(df_virhost.xref_B)),
                            ['xref_A', 'xref_B']]).size)

