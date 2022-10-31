#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

from sklearn import svm
from sklearn.decomposition import PCA
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.preprocessing import MinMaxScaler

from matplotlib import pyplot
import matplotlib.pyplot as plt
import plotly.express as px
import xgboost as xgb


# In[2]:


countsRaw = pd.read_csv('/Users/elikond/Downloads/seq_counts.cts', sep='\t')
engraftmentRaw = pd.read_csv('/Users/elikond/Downloads/Engraftment_Sheet.tsv', sep='\t')


# In[3]:


countsRaw.dropna(subset=['Genes'], inplace = True)


# In[4]:


unlabelled_patients = ['14415-30.', '14415-04', '14415-26']
cohort_2 = ['14415-22', '14415-12', '14415-13', '14415-14', '14415-16']

all_cohort_counts = countsRaw.drop(unlabelled_patients, axis=1)
two_cohort_counts = countsRaw.drop(unlabelled_patients + cohort_2, axis=1)
all_cohort_counts.head()


# In[5]:


def process_engraftment(engraftment_df):
    new_engraftment = engraftment_df.copy()
    new_engraftment['Lowercase_Engraftment'] = new_engraftment.Engraftment.str.lower()
    new_engraftment['Binarized Engraftment'] = new_engraftment.Lowercase_Engraftment.eq('h').mul(1)
    new_engraftment.set_index('Patient', drop = True, inplace = True)
    
    return new_engraftment

processed_engraftment = process_engraftment(engraftmentRaw)

all_cohort_engraftment = processed_engraftment[~processed_engraftment.index.isin(unlabelled_patients)]
two_cohort_engraftment = processed_engraftment[~processed_engraftment.index.isin(unlabelled_patients + cohort_2)]
all_cohort_engraftment.head()


# In[6]:


#Printing ratio of low engraftment in engraftment dataframe
def find_percent_zero(engraftment_df, type):
    num_zeroes = engraftment_df['Binarized Engraftment'][engraftment_df['Binarized Engraftment'] == 0].count()
    percent_zeroes = num_zeroes/len(engraftment_df)
    print(f'Ratio of {type} Patients with Low Engraftment: %0.3f' %percent_zeroes)

find_percent_zero(all_cohort_engraftment, 'All')
find_percent_zero(two_cohort_engraftment, 'Cohort 1 and 3')


# In[7]:


def mean_norm(df_input):
    return df_input.apply(lambda x: (x-x.mean())/ x.std(), axis=0)


# In[8]:


#Standardizing by rows and transposing
def process_counts(input_df):
    temp_df = mean_norm(input_df.iloc[:,1:])
    tranposed_data = temp_df.T
    tranposed_data.columns = list(input_df['Genes'])

    return tranposed_data

ProcessedCounts_all_cohort = process_counts(all_cohort_counts)
ProcessedCounts_two_cohort = process_counts(two_cohort_counts)

ProcessedCounts_all_cohort.head()


# In[9]:


#Dimensionality reduction and adding in engraftment data

def dim_reduction(prepped_data, engraftment_df):
    pca = PCA(n_components=len(prepped_data))
    reduced_data = pca.fit_transform(prepped_data)
    components = pca.components_
    pca_df = pd.DataFrame(reduced_data)
    pca_df.index = prepped_data.index
    pca_df = pca_df.add_prefix('PCA_')
    final_df = pca_df.join(engraftment_df['Binarized Engraftment'])
    return pca, final_df

all_cohort_pca, FinalData_all_cohort = dim_reduction(ProcessedCounts_all_cohort, all_cohort_engraftment)
two_cohort_pca, FinalData_two_cohort = dim_reduction(ProcessedCounts_two_cohort, two_cohort_engraftment)
FinalData_all_cohort.head()

