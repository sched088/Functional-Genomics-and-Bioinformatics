import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import cluster


# Set variables
fname = 'Data1.npz'
k_list = [10,20,50]

# Load data
data_1 = np.load(fname)

# Understand the data
print(data_1.files)
print('Length of Gene_Name list: ' + str(len(data_1['Gene_Name'])))
print('Length of P_value list: ' + str(len(data_1['P_value'])))
print('Length of SeqData list: ' + str(len(data_1['SeqData'])))
    # Confirmed length is the same for each file

# plt.hist(data_1['SeqData'], bins = 20)
# plt.show()
# plt.hist(data_1['SeqData'][:1000], bins = 20)
# plt.show()


# Convert npz to dataframe for future data manipulation
d1_df = pd.DataFrame.from_dict({item: data_1[item] for item in data_1.files}, orient = 'index')
d1_df = d1_df.T

# print(d1_df['SeqData'].iloc[:10])

# functional enrichment analysis
# Extract first 200 genes from Gene_Name list:
def fea(data_1): 
    np.savetxt('200_Gene_Names.txt', [data_1['Gene_Name'][:200]], fmt='%s', delimiter=',')


# Extract first 1000 genes from Gene_Name list for kmeans:
# k=10, 20, 50
# For each case, plot the histogram of the cluster sizes.
def kmeans(seqdata, genelabels, ks):
    # object_count = 1000
    # cluster_on = 'SeqData'
    # fit_data = np.asmatrix(data[cluster_on][:object_count])
    # print(np.asarray(fit_data))
    for k in ks:
        d1_kmeans = cluster.KMeans(n_clusters=k, random_state=0).fit(seqdata)

        if k == 20:
            kmeans_df = pd.DataFrame()
            kmeans_df['cluster_idx'] = d1_kmeans.labels_
            kmeans_df['data_values'] = genelabels

        # plt.hist(d1_kmeans.labels_, bins = k)
        # plt.xticks(np.arange(k))
        # plt.show()

    return kmeans_df

def kmeans_fea(kmeans_data):
    selected_cluster = 12
    np.savetxt('kmeans_fea.txt', [kmeans_data.loc[kmeans_data['cluster_idx'] == selected_cluster, 'data_values']], fmt='%s', delimiter=',')

    

kmeans_df = kmeans(data_1['SeqData'][:1000], data_1['Gene_Name'][:1000], k_list)
kmeans_fea(kmeans_df)