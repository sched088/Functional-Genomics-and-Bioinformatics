'''
p <= i/N * (alpha)
alpha = p*N/i

Assumptions:
 p <= 0.05
 N = 2239 (from problem 2)
'''

from statsmodels.stats import multitest
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import sklearn
from sklearn.feature_selection import SelectFdr

load_path_HiSeqV2 = "/Users/sched/Downloads/BIOC HW2/SeqData.txt"
load_path_HT_HG = "/Users/sched/Downloads/BIOC HW2/ArrayData.txt"


df_seqdata = pd.read_csv(load_path_HiSeqV2, sep='\t')
df_arraydata = pd.read_csv(load_path_HT_HG, sep='\t')


df_combined = pd.merge(df_seqdata, df_arraydata, on=['Sample ID', 'Group Label'], how='outer')

# Filter genes with 0 expression across all patients
df = df_combined.loc[:, (df_combined != 0).any(axis=0)]

print(df_combined.shape)
print(df.shape)



# Split into Groups
df_group1 = df[df['Group Label'] == "Group 1"]
df_group2 = df[df['Group Label'] == "Group 2"]

df_group1 = df_group1.head(10)
df_group2 = df_group2.head(10)

ttest_results = [[header] for header in df.iloc[:,2:].columns]
# idx = df_group1.index.intersection(df_group2.index)
# print(ttest_ind(df_group1.iloc[:,2:], df_group2.iloc[:,2:]))
ttest, pvalue = ttest_ind(df_group1.iloc[:,2:], df_group2.iloc[:,2:])

filter = pvalue < 0.05
filtered_pvalues = pvalue[filter]

n=0
for _ in ttest_results:
    _.append(ttest[n])
    _.append(pvalue[n])
    n+1


bon_alpha = 0.05 / len(filtered_pvalues)
reject, pvals_corrected, alphaSidak, alphacBong = multitest.multipletests(filtered_pvalues, alpha=bon_alpha, method='bonferroni', is_sorted=False, returnsorted=False)
# print(reject[:10])
print("Bonferroni # rejects: " + str(sum(reject)))
# print((pvals_corrected))

# # i = 20 // a = 5.6 // 0.05 * 20 / 2239 = 0.000447
# i = 20
# fdr_alpha = (0.05 * i) / len(filtered_pvalues)
# rejected, pvalue_corrected = multitest.fdrcorrection(filtered_pvalues,alpha=fdr_alpha)
# print("FDR 20 # rejects: " + str(sum(rejected)))

# # i = 50 // a = 2.24 // 0.00112
# i = 50
# fdr_alpha = (0.05 * i) / len(filtered_pvalues)
# rejected, pvalue_corrected = multitest.fdrcorrection(filtered_pvalues,alpha=fdr_alpha)
# print("FDR 50 # rejects: " + str(sum(rejected)))

# # i = 100 // a = 1.12 // 0.0022
# i = 100
# fdr_alpha = (0.05 * i) / len(filtered_pvalues)
# rejected, pvalue_corrected = multitest.fdrcorrection(filtered_pvalues,alpha=fdr_alpha)
# print("FDR 100 # rejects: " + str(sum(rejected)))

# # i = 200 // a = 0.56 // 0.0044
# i = 200
# fdr_alpha = (0.05 * i) / 300     # len(filtered_pvalues)
# rejected, pvalue_corrected = multitest.fdrcorrection(filtered_pvalues,alpha=fdr_alpha)
# print("FDR 200 # rejects: " + str(sum(rejected)))
