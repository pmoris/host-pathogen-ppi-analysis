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

