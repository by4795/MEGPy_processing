import seaborn as sns
import numpy as np
import csv
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy import stats


data = pd.read_csv("/home/ir-sidd1/rds/rds-ukaea-ap002-mOlK9qn0PlQ/ir-sidd1/TGLF/tglf_nn_gyrokit/outputs/tglf_inputs_aaron.csv") #dataframe

chosen_params = ["KAPPA_LOC", "S_KAPPA_LOC", "DRMAJDX_LOC", "DELTA_LOC", "S_DELTA_LOC", "RMIN_LOC"] #geometric parameters we care about
#should already be preselected however

print(data.size)

#now what i want to do is to go through this data frame and remove the most extreme values
#i could either replace these with the mean value or remove the row of the data

acceptable_std = 1

for parameter in chosen_params:
    #find mean
    #q = data[parameter].quantile(0.99)
    #data[data[parameter] < q]
    #find 99th percentile
    #data[(np.abs(stats.zscore(data)) < acceptable_std ).all(axis=1)]
    q_low = data[parameter].quantile(0.01)
    q_hi  = data[parameter].quantile(0.99)

    data = data[(data[parameter] < q_hi) & (data[parameter] > q_low)]

print(data.size)

fig = sns.pairplot(data, corner=True)

fig.savefig("corner_plot_geometries_outliers_removed_aaron.png")

