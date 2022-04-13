from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set variables
fname = 'Data2.npz'

# Load data
data_2 = np.load(fname)

# Understand the data
# Uncomment for data overview
# print(data_2.files)
# for file in data_2.files:
#     print('length of {} list: {}'.format(str(file), str(len(data_2[file]))))
#     print(data_2[file][:10])


# Run KNN classifier with k=1, 3, 5 
# Use the samples in training_data as training data and the samples in testing_data as test data. 
# Report your classification accuracy for the three cases.
def knn_1(data, ks):
    for k in ks:
        knn_k = KNeighborsClassifier(k)
        knn_fit = knn_k.fit(data['training_data'], data['training_label'])
        knn_pred = knn_k.predict(data['testing_data'])
        knn_score = knn_k.score(data['testing_data'], data['testing_label'])

        # print(knn_k)
        # print(knn_fit)
        # print(knn_pred)
        # print(knn_score)
        print('Accuracy of testing_data using training_data and k={} is {}%'.format(str(k), str(knn_score*100)))

def knn_2(data, ks, flag):
        if flag == 1:
            scaler = StandardScaler()
            scaler.fit(data['training_data']) 
            for k in ks:
                knn_k = KNeighborsClassifier(k)
                knn_k.fit(scaler.transform(data['training_data']), data['training_label'])
                knn_score = knn_k.score(scaler.transform(data['testing_data']), data['testing_label'])
                print('Accuracy of testing_data using {} training_data and k={} is {}%'.format('all', str(k), str(format(knn_score*100,".2f"))))

        if flag == 0:
            scaler = StandardScaler()
            scaler.fit(data['training_data'][:,:1000]) 
            for k in ks:
                knn_k = KNeighborsClassifier(k)
                knn_k.fit(scaler.transform(data['training_data'][:,:1000]), data['training_label'])
                knn_score = knn_k.score(scaler.transform(data['testing_data'][:,:1000]), data['testing_label'])
                print('Accuracy of testing_data using {} training_data and k={} is {}%'.format(' the first 1000 points of', str(k), str(format(knn_score*100,".2f"))))

        else:
            return "please set flag to 1 or 0 only"
        # print(knn_k)
        # print(knn_fit)
        # print(knn_pred)
        # print(knn_score)

k_list = [1,3,5]
flag = 1
knn_2(data_2, k_list, flag)
flag = 0
knn_2(data_2, k_list, flag)