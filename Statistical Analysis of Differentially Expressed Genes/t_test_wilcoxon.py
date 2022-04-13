import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, wilcoxon
import matplotlib.pyplot as plt

load_path_HiSeqV2 = "/Users/sched/Downloads/BIOC HW2/SeqData.txt"
load_path_HT_HG = "/Users/sched/Downloads/BIOC HW2/ArrayData.txt"


df_seqdata = pd.read_csv(load_path_HiSeqV2, sep='\t')
df_arraydata = pd.read_csv(load_path_HT_HG, sep='\t')

# new_index = df_seqdata["Sample ID"] 
# df_seqdata.index = new_index

# new_index = df_arraydata["Sample ID"] 
# df_arraydata.index = new_index

# print(df_seqdata.head(10))
# print(df_seqdata.tail(10))

# print(df_arraydata.head(10))
# print(df_arraydata.tail(10))


# df_combined = pd.concat([df_seqdata, df_arraydata], axis=0)

df_combined = pd.merge(df_seqdata, df_arraydata, on=['Sample ID', 'Group Label'], how='outer')

# Filter genes with 0 expression across all patients
df = df_combined.loc[:, (df_combined != 0).any(axis=0)]


# print(df_combined.head(10))
# print(df_combined.tail(10))

# print(df_seqdata.shape)
# print(df_arraydata.shape)
print(df_combined.shape)
print(df.shape)


# get union of data
# df_combined.dropna(axis = 1)

# Split into Groups
df_group1 = df[df['Group Label'] == "Group 1"]
df_group2 = df[df['Group Label'] == "Group 2"]

df_group1 = df_group1.head(10)
df_group2 = df_group2.head(10)

ttest_results = [[header] for header in df.iloc[:,2:].columns]
# idx = df_group1.index.intersection(df_group2.index)
# print(ttest_ind(df_group1.iloc[:,2:], df_group2.iloc[:,2:]))
ttest, pvalue = ttest_ind(df_group1.iloc[:,2:], df_group2.iloc[:,2:])

n=0
for _ in ttest_results:
    _.append(ttest[n])
    _.append(pvalue[n])
    n+=1

sorted_results = sorted(ttest_results, key=lambda x: x[2])

filter = pvalue < 0.05
filtered_pvalues = pvalue[filter]

print("T-Test Results - The 10 genes with the most significant p-values are as follows [gene, t-stat, p-value]: " + str(sorted_results[:10]))
print("T-Test Results - The 10 genes with the most significant p-values are as follows [gene, t-stat, p-value]: " + str(sorted_results[-10:]))


print("There are " + str(len(filtered_pvalues)) + " genes with significant p-values")

# w_ttest, w_pvalue = wilcoxon(df_group1.iloc[:,2:], df_group2.iloc[:,2:])

wilcox_results = [[header] for header in df.iloc[:,2:].columns]
w_pvalue_array = []

n=0
for result in wilcox_results:
    w_ttest, w_pvalue = wilcoxon(df_group1.iloc[:,n+2], df_group2.iloc[:,n+2])
    result.append(w_ttest)
    result.append(w_pvalue)
    w_pvalue_array.append(w_pvalue)
    n+=1

sorted(wilcox_results, key=lambda x: x[2])

print("Wilcoxon Results - The genes with the top 10 p-values are as follows [gene, t-stat, p-value]: " + str(wilcox_results[-10:]))

p_bp_ttest = plt.hist(filtered_pvalues, label= "ttest",
                     density= True,
                     alpha=0.75)
plt.show()
p_bp_wilcox = plt.hist(w_pvalue_array, label= "wilcox",
                       density= True,
                       alpha=0.75)
plt.show()
