# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:22:01 2016

@author: bmanubay
"""


import numpy as np
import pandas as pd

# Read in property data produced from puresrch.py and binsrch.py
pathdfp = "/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/"
pathdfm = "/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/"

a1 = pd.read_csv(pathdfp+"dens_pure.csv", sep=';')   
a2 = pd.read_csv(pathdfp+"hmol_pure.csv", sep=';')
a3 = pd.read_csv(pathdfp+"cpmol_pure.csv", sep=';')

b1 = pd.read_csv(pathdfm+"dens_bin.csv", sep=';')
b2 = pd.read_csv(pathdfm+"eme_bin.csv", sep=';')
b3 = pd.read_csv(pathdfm+"emcp_bin.csv", sep=';')

# Create new columns for combined SMILES strings on mixtures
b1['SMILES_mix'] = b1['SMILES1'].map(str) + '__' + b1['SMILES2'].map(str)
b2['SMILES_mix'] = b2['SMILES1'].map(str) + '__' + b2['SMILES2'].map(str)
b3['SMILES_mix'] = b3['SMILES1'].map(str) + '__' + b3['SMILES2'].map(str)

# Merge property sets using filename as index
c1 = pd.merge(a1,b1,how='outer',on=['filename'])
c2 = pd.merge(a2,b2,how='outer',on=['filename'])
c3 = pd.merge(a3,b3,how='outer',on=['filename'])

# Find unique compoounds/mixtures for each journal in merged datasets 
d1 = c1.groupby('filename').SMILES.unique()
d1 = d1.reset_index()
d1.columns = ['filename', 'Pure_Dens_Molecules']
e1 = c1.groupby('filename').SMILES_mix.unique()
e1 = e1.reset_index()
e1.columns = ['filename', 'Bin_Dens_Mixtures']
d2 = c2.groupby('filename').SMILES.unique()
d2 = d2.reset_index()
d2.columns = ['filename', 'Pure_Hmol_Molecules']
e2 = c2.groupby('filename').SMILES_mix.unique()
e2 = e2.reset_index()
e2.columns = ['filename', 'Bin_EME_Mixtures']
d3 = c3.groupby('filename').SMILES.unique()
d3 = d3.reset_index()
d3.columns = ['filename', 'Pure_Cpmol_Molecules']
e3 = c3.groupby('filename').SMILES_mix.unique()
e3 = e3.reset_index()
e3.columns = ['filename', 'Bin_EMCp_Mixtures']

# Merge the unique component/mixture list with filename as index
f1 = pd.merge(d1,e1,how='outer',on=['filename'])
f2 = pd.merge(d2,e2,how='outer',on=['filename'])
f3 = pd.merge(d3,e3,how='outer',on=['filename'])

g1 = pd.merge(f1,f2,how='outer',on=['filename'])
g2 = pd.merge(g1,f3,how='outer',on=['filename'])

# Convert data type to string to allow simpler application of counting process
g2 = g2.astype('str')

# Format strings in columns
g2['Pure_Dens_Molecules'] = g2['Pure_Dens_Molecules'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Bin_Dens_Mixtures'] = g2['Bin_Dens_Mixtures'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Pure_Hmol_Molecules'] = g2['Pure_Hmol_Molecules'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Bin_EME_Mixtures'] = g2['Bin_EME_Mixtures'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Pure_Cpmol_Molecules'] = g2['Pure_Cpmol_Molecules'].map(lambda x: x.lstrip('"[').rstrip(']"'))
g2['Bin_EMCp_Mixtures'] = g2['Bin_EMCp_Mixtures'].map(lambda x: x.lstrip('"[').rstrip(']"'))

# Replace string nan with np.nan float
g2.replace('nan', np.nan, inplace=True)
g2.replace('[nan]', np.nan, inplace=True)

cols = ['Pure_Dens_Molecules','Bin_Dens_Mixtures','Pure_Hmol_Molecules','Bin_EME_Mixtures','Pure_Cpmol_Molecules','Bin_EMCp_Mixtures']

# Return boolean data frame to count how many properties are measured per journal
g3 = pd.isnull(g2)
g2['Property_Score'] = g3[cols].apply(lambda x: (x==False).sum(), axis=1)

# Count ocurrences of nan per row in order to be able to subtract from component counts later
g2['nan_Count'] = g3[cols].apply(lambda x: (x==True).sum(), axis=1)

# Turn the dataframe back to string format and keep the count columns as float
g2 = g2.astype('str')
g2['Property_Score'] = g2.Property_Score.astype('float')
g2['nan_Count'] = g2.nan_Count.astype('float')
g2['nan_Count'] = -g2.nan_Count

# Count number of components/mixtures per property in all journals (counts nan as well, hence why the previous nan count is saved)
g2['Pure_Dens_Molecules'] = g2['Pure_Dens_Molecules'].str.split()
g2['Pure_Dens_Counts'] = g2['Pure_Dens_Molecules'].map(lambda x: len(x))
g2['Bin_Dens_Mixtures'] = g2['Bin_Dens_Mixtures'].str.split()
g2['Bin_Dens_Counts'] = g2['Bin_Dens_Mixtures'].map(lambda x: len(x))
g2['Pure_Hmol_Molecules'] = g2['Pure_Hmol_Molecules'].str.split()
g2['Pure_Hmol_Counts'] = g2['Pure_Hmol_Molecules'].map(lambda x: len(x))
g2['Bin_EME_Mixtures'] = g2['Bin_EME_Mixtures'].str.split()
g2['Bin_EME_Counts'] = g2['Bin_EME_Mixtures'].map(lambda x: len(x))
g2['Pure_Cpmol_Molecules'] = g2['Pure_Cpmol_Molecules'].str.split()
g2['Pure_Cpmol_Counts'] = g2['Pure_Cpmol_Molecules'].map(lambda x: len(x))
g2['Bin_EMCp_Mixtures'] = g2['Bin_EMCp_Mixtures'].str.split()
g2['Bin_EMCp_Counts'] = g2['Bin_EMCp_Mixtures'].map(lambda x: len(x))


cols = ['Pure_Dens_Counts','Bin_Dens_Counts','Pure_Hmol_Counts','Bin_EME_Counts','Pure_Cpmol_Counts','Bin_EMCp_Counts', 'nan_Count']

# Create column which is count of number of components/mixtures per journal
g2['Component_Score'] = g2[cols].sum(axis=1)

cols = ['Component_Score', 'Property_Score']

# Create aggregate score based on the property and component scores 
g2['Aggregate_Score'] = g2[cols].sum(axis=1)

# Sort data frame by aggregate score in descending order 
g2 = g2.sort_values(by='Aggregate_Score', ascending=False)

# Save output data frame as csv and pickle
g2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Journal_scores_for_Ken.csv", sep =';')
g2.to_pickle("/home/bmanubay/.thermoml/tables/Ken/Journal_scores_for_Ken.pkl")