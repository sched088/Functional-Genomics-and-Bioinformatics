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


def svc_analysis(data):
    scaler = StandardScaler()
    scaler.fit(data['training_data'][:,:1000]) 
    svc = LinearSVC(random_state=0)
    svc.fit(scaler.transform(data['training_data'][:,:1000]), data['training_label'])
    svc_score = svc.score(scaler.transform(data['testing_data'][:,:1000]), data['testing_label'])

    print('Accuracy of testing_data using {} training_data and SVC approach is {}%'.format(' the first 1000 points of', str(format(svc_score*100,".2f"))))

svc_analysis(data_2)
