#Importing all necessary packages
import pandas as pd
from sklearn.decomposition import PCA
import plotly.express as px
from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import StratifiedGroupKFold
import numpy as np
from numpy import mean
from numpy import isnan
from numpy import asarray
from numpy import polyfit
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
from sklearn.model_selection import LeaveOneOut

from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import ExtraTreeClassifier
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

from plotly.offline.offline import plot

#Reading in count + engraftment data
count_df = pd.read_csv('../counts.cts', sep='\t')
unprocessed_engraftment_df = pd.read_csv('../Engraftment_Sheet2.tsv', sep='\t')

def prep_engraftment_data(engraftment_data):
    engraftment_data['Lowercase_Engraftment'] = engraftment_data.Engraftment.str.lower()
    engraftment_data['Binarized Engraftment'] = engraftment_data.Lowercase_Engraftment.eq('h').mul(1)
    engraftment_data.set_index('Patient', drop = True, inplace = True)

    num_zeroes = engraftment_data['Binarized Engraftment'][engraftment_data['Binarized Engraftment'] == 0].count()
    percent_zeroes = num_zeroes/len(engraftment_data)
    print('Ratio of patients with low engraftment: %0.3f' %percent_zeroes)

    return engraftment_data

engraftment_df = prep_engraftment_data(unprocessed_engraftment_df)

def mean_norm(df_input):
    return df_input.apply(lambda x: (x-x.mean())/ x.std(), axis=0)


def transpose(df, count_data = count_df):
  if 'Genes' in df.columns:
    df.drop('Genes', axis = 1, inplace = True)
  tranposed_data = df.T
  tranposed_data.columns = list(count_data['Genes'])
  return tranposed_data

def prep_count_data(count_data):
    #Removing patients for which we do not have engraftment labels
    labeled_count_df = count_data.drop(['14415-30.', '14415-04', '14415-26'], axis=1)
    standardized_df = mean_norm(labeled_count_df.iloc[:,1:])
    prepped_data = transpose(standardized_df)
    pca = PCA(n_components=24)
    reduced_data = pca.fit_transform(prepped_data)

    pca_df = pd.DataFrame(reduced_data)
    pca_df.index = prepped_data.index
    pca_df = pca_df.add_prefix('PCA_')

    fig = px.imshow(pca_df, labels=dict(x="Dimensionally Reduced Genes", y="Patients"))
    fig.show()

    pca_components_genes = pd.DataFrame(pca.components_,columns=prepped_data.columns,index = pca_df.columns)

    return labeled_count_df, pca_df, pca_components_genes

labeled_count_df, reduced_data, pca_df, pca_components_genes = prep_count_data(count_df)

#i is the pca component index
def individual_pca_heatmap(count_data, pca_gene_arr, lim, i):
  #if out limit it less than one, we are defining a threshold for the genes to pass
  if lim < 1:
    s = pca_gene_arr.iloc[i,:]
    top_genes = list(s[s > lim].index)
  #if our limit is greater than one, we are defining the number of genes to chose
  if lim > 1:
    top_genes = list(pca_gene_arr.iloc[i,:].nlargest(lim).index)
  top_genes_counts = count_data[count_data.Genes.isin(top_genes)]
  top_genes_counts.set_index('Genes', drop = True, inplace = True)
  plot_data = mean_norm(top_genes_counts)
  #print(plot_data)
  fig = px.imshow(plot_data, labels=dict(x="Patients", y="Genes"), title = 'PCA Component ' + str(i))
  fig.show()

individual_pca_heatmap(labeled_count_df, pca_components_genes, 15, 0)

def run_svm(pca_df, engraftment_df, reduced_data):
    #Combining PCA and Engraftment dataframes to get the patients in the same order
    combined_pca_engraftment = pd.concat([pca_df, engraftment_df], axis = 1)
    engraftment_labels = combined_pca_engraftment['Binarized Engraftment']

    #Checking that the order of patients was not changed from the reduced_data after dfs were combined
    print(combined_pca_engraftment.index == pca_df.index)

    clf = svm.SVC(kernel='rbf', C=1, random_state=42)
    cv = StratifiedGroupKFold(n_splits=3)
    #group = combined_pca_engraftment['Binarized Engraftment'].astype('str') + '_' +  combined_pca_engraftment['Cohort'].astype('str')
    group = combined_pca_engraftment['Cohort']
    scores = cross_val_score(clf, reduced_data, engraftment_labels, groups = group, cv=cv, scoring='accuracy')
    print(scores)
    print("%0.2f accuracy with a standard deviation of %0.2f" % (scores.mean(), scores.std()))

    for train_index, test_index in cv.split(reduced_data, engraftment_labels, group):
        print("TRAIN:", train_index, "TEST:", test_index)
