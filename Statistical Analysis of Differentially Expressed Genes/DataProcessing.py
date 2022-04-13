"""
Preprocessing Spreadsheets
"""

import pandas as pd
import numpy as np
import random 
import matplotlib.pyplot as plt


load_path_HiSeqV2 = "/Users/sched/Downloads/BIOC HW2/HiSeqV2"
load_path_HT_HG = "/Users/sched/Downloads/BIOC HW2/HT_HG-U133A"
save_path_HiSeqV2 = "/Users/sched/Downloads/BIOC HW2/SeqData.txt"
save_path_HT_HG = "/Users/sched/Downloads/BIOC HW2/ArrayData.txt"

def preprocessing(load_path_HiSeqV2, load_path_HT_HG, save_path_HiSeqV2, save_path_HT_HG):
    df_hiseq = pd.read_csv(load_path_HiSeqV2, sep='\t')
    df_hiseq = df_hiseq.T
    new_header = df_hiseq.iloc[0] 
    df_hiseq = df_hiseq[1:] 
    df_hiseq.columns = new_header 
    df_hiseq.reset_index(inplace=True)
    df_hiseq = df_hiseq.rename(columns = {'index':'Sample ID'})

    df_hthg = pd.read_csv(load_path_HT_HG, sep='\t')
    df_hthg = df_hthg.T
    new_header = df_hthg.iloc[0] 
    df_hthg = df_hthg[1:] 
    df_hthg.columns = new_header 
    df_hthg.reset_index(inplace=True)
    df_hthg = df_hthg.rename(columns = {'index':'Sample ID'})

    array_hthg = np.array(df_hthg.columns)
    array_hiseq = np.array(df_hiseq.columns)
    gene_list = list(set(list(array_hiseq) + list(array_hthg)))

    df_clinical = pd.read_csv("ov_tcga_clinical_data.tsv", sep='\t')

    # Group the data by deceased at 36 months
    def group_label(row):
        if row["Overall Survival (Months)"] < 36 and row["Overall Survival Status"] == "1:DECEASED":
            return "Group 1"
        if row["Overall Survival (Months)"] > 36 and row["Overall Survival Status"] == "0:LIVING":
            return "Group 2"
        else:
            return "Group 0"

    df_clinical['Group Label'] = df_clinical.apply (lambda row: group_label(row), axis=1)

    df_clin_hthg = pd.merge(df_clinical[['Group Label', 'Sample ID']], df_hthg, on=['Sample ID'])
    df_clin_hiseq = pd.merge(df_clinical[['Group Label', 'Sample ID']], df_hiseq, on=['Sample ID'])

    df_clin_hthg = df_clin_hthg[df_clin_hthg['Group Label'] != 'Group 0']
    df_clin_hiseq = df_clin_hiseq[df_clin_hiseq['Group Label'] != 'Group 0']

    # print(df_clin_hiseq.shape)
    # print(df_clin_hthg.shape)

    df_inner_intersect = pd.concat([df_clin_hiseq, df_clin_hthg], axis=0, join='inner')
    df_clin_hthg = pd.concat([df_inner_intersect, df_clin_hthg], axis=0, join='inner')
    df_clin_hiseq = pd.concat([df_inner_intersect, df_clin_hiseq], axis=0, join='inner')
    df_inner_intersect.dropna(axis=1)
    df_clin_hthg.dropna(axis=1)
    df_clin_hiseq.dropna(axis=1)

    df_clin_hthg.to_csv(save_path_HT_HG, index_label=False, sep ='\t')
    df_clin_hiseq.to_csv(save_path_HiSeqV2, index_label=False, sep ='\t')

    # print(df_inner_intersect.shape)
    # print(df_clin_hiseq.shape)
    # print(df_clin_hthg.shape)

    return df_inner_intersect, gene_list

def plot_overlap(df_inner_intersect, gene_list):
    gene_count = 0
    points_to_plot = []
    while gene_count <= 1000:
        gene_samples = random.sample(gene_list, gene_count)	
        shared_gene_count = len(list(set(df_inner_intersect).intersection(gene_samples)))
        coordinates = (len(gene_samples), shared_gene_count)
        points_to_plot.append(coordinates)
        gene_count += 10

    zip(*points_to_plot)
    plt.scatter(*zip(*points_to_plot))
    plt.show()

df_inner_intersect, gene_list = preprocessing(load_path_HiSeqV2, load_path_HT_HG, save_path_HiSeqV2, save_path_HT_HG)
plot_overlap(df_inner_intersect, gene_list)


