# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:09:31 2016

@author: bmanubay
"""

import cirpy
import numpy as np
import pandas as pd
from sklearn.externals.joblib import Memory

mem = Memory(cachedir="/home/bmanubay/.thermoml/")

@mem.cache
def resolve_cached(x, rtype):
   return cirpy.resolve(x, rtype)

# Define list of all alkane SMILES strings that appear in all of our data
SMILES_alk = ['C', 'CC', 'CCC', 'CCCC', 'CCCCC', 'CCCCCC', 'CCCCCCC', 'CCCCCCCC', 'CCCCCCCCC', 'CCCCCCCCCC', 'CC(C)C', 'CCC(C)C', 'CCCC(C)C', 'C1CCCCC1', 'CC1CCCCC1', 'CCCCCC(C)C', 'CC(C)C(C)C', 'CCC(C)(C)C', 'CCC(C)CC', 'CCCC(C)C', 'CC(C)CC(C)(C)C', 'C1CCCC1', 'C1CCCCCCC1', 'CCC1CCCCC1', 'CC1CCC(C)C(C)C1', 'CCCCC1CCCCC1', 'CC1CCCC1', 'CCCCCCC(C)C', 'CCCCCCCC(C)C', 'CCCCC(C)CCC', 'CCC(C)CCC(C)CC', 'CCC(C)CC(C)CC', 'CCCCCC(C)C', 'C1CCCCCC1', 'CC(C)C(C)C', 'CCC(C)(C)C']

S = pd.DataFrame({'SMILES': SMILES_alk}, columns = ['SMILES'])

# Binary mixtures
aa1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/actcoeff_bin.csv", sep=';', index_col= 'Unnamed: 0')
aa2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/dens_bin.csv", sep=';', index_col= 'Unnamed: 0')
aa3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/dielec_bin.csv", sep=';', index_col= 'Unnamed: 0')
aa4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/eme_bin.csv", sep=';', index_col= 'Unnamed: 0')
aa5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/emcp_bin.csv", sep=';', index_col= 'Unnamed: 0')
aa6 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/emv_bin.csv", sep=';', index_col= 'Unnamed: 0')
aa7 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/sos_bin.csv", sep=';', index_col= 'Unnamed: 0')


# Binary Mixtures with alkane-alkane mixtures removed
cc1c = pd.concat([aa1["x1"], aa1["x2"], aa1['SMILES1'], aa1['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"]) 
ind = np.logical_not(cc1c.SMILES1.isin(S.SMILES) & cc1c.SMILES2.isin(S.SMILES))
cc1c = cc1c[ind]
cc1c = cc1c.drop(['x1','x2'], axis=1)
count1c = pd.Series(cc1c.squeeze().values.ravel()).value_counts()
count1c = count1c.reset_index()
count1c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc1 = pd.concat([count1c["SMILES"], count1c["Count"]], axis=1, keys=["SMILES", "Count"])
cc2c = pd.concat([aa2["x1"], aa2["x2"], aa2['SMILES1'], aa2['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"]) 
ind = np.logical_not(cc2c.SMILES1.isin(S.SMILES) & cc2c.SMILES2.isin(S.SMILES))
cc2c = cc2c[ind]
cc2c = cc2c.drop(['x1','x2'], axis=1)
count2c = pd.Series(cc2c.squeeze().values.ravel()).value_counts()
count2c = count2c.reset_index()
count2c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc2 = pd.concat([count2c["SMILES"], count2c["Count"]], axis=1, keys=["SMILES", "Count"])
cc3c = pd.concat([aa3["x1"], aa3["x2"], aa3['SMILES1'], aa3['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"]) 
ind = np.logical_not(cc3c.SMILES1.isin(S.SMILES) & cc3c.SMILES2.isin(S.SMILES))
cc3c = cc3c[ind]
cc3c = cc3c.drop(['x1','x2'], axis=1)
count3c = pd.Series(cc3c.squeeze().values.ravel()).value_counts()
count3c = count3c.reset_index()
count3c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc3 = pd.concat([count3c["SMILES"], count3c["Count"]], axis=1, keys=["SMILES", "Count"])
cc4c = pd.concat([aa4["x1"], aa4["x2"], aa4['SMILES1'], aa4['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"]) 
ind = np.logical_not(cc4c.SMILES1.isin(S.SMILES) & cc4c.SMILES2.isin(S.SMILES))
cc4c = cc4c[ind]
cc4c = cc4c.drop(['x1','x2'], axis=1)
count4c = pd.Series(cc4c.squeeze().values.ravel()).value_counts()
count4c = count4c.reset_index()
count4c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc4 = pd.concat([count4c["SMILES"], count4c["Count"]], axis=1, keys=["SMILES", "Count"])
cc5c = pd.concat([aa5["x1"], aa5["x2"], aa5['SMILES1'], aa5['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"])
ind = np.logical_not(cc5c.SMILES1.isin(S.SMILES) & cc5c.SMILES2.isin(S.SMILES))
cc5c = cc5c[ind] 
cc5c = cc5c.drop(['x1','x2'], axis=1)
count5c = pd.Series(cc5c.squeeze().values.ravel()).value_counts()
count5c = count5c.reset_index()
count5c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc5 = pd.concat([count5c["SMILES"], count5c["Count"]], axis=1, keys=["SMILES", "Count"])
cc6c = pd.concat([aa6["x1"], aa6["x2"], aa6['SMILES1'], aa6['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"]) 
ind = np.logical_not(cc6c.SMILES1.isin(S.SMILES) & cc6c.SMILES2.isin(S.SMILES))
cc6c = cc6c[ind]
cc6c = cc6c.drop(['x1','x2'], axis=1)
count6c = pd.Series(cc6c.squeeze().values.ravel()).value_counts()
count6c = count6c.reset_index()
count6c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc6 = pd.concat([count6c["SMILES"], count6c["Count"]], axis=1, keys=["SMILES", "Count"])
cc7c = pd.concat([aa7["x1"], aa7["x2"], aa7['SMILES1'], aa7['SMILES2']], axis=1, keys = ["x1", "x2", "SMILES1", "SMILES2"]) 
ind = np.logical_not(cc7c.SMILES1.isin(S.SMILES) & cc7c.SMILES2.isin(S.SMILES))
cc7c = cc7c[ind]
cc7c = cc7c.drop(['x1','x2'], axis=1)
count7c = pd.Series(cc7c.squeeze().values.ravel()).value_counts()
count7c = count7c.reset_index()
count7c.rename(columns={"index":"SMILES",0:"Count"},inplace=True)
ccc7 = pd.concat([count7c["SMILES"], count7c["Count"]], axis=1, keys=["SMILES", "Count"])


# All data counts
c1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_cpmol.csv", sep=';', usecols=['SMILES', 'Count'])
c2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_dens.csv", sep=';', usecols=['SMILES', 'Count'])
c3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_dielec.csv", sep=';', usecols=['SMILES', 'Count'])
c4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_hmol.csv", sep=';', usecols=['SMILES', 'Count'])
c5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_sos.csv", sep=';', usecols=['SMILES', 'Count'])

cc1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_actcoeff.csv", sep=';', usecols=['SMILES', 'Count'])
cc2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_dens.csv", sep=';', usecols=['SMILES', 'Count'])
cc3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_dielec.csv", sep=';', usecols=['SMILES', 'Count'])
cc4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_eme.csv", sep=';', usecols=['SMILES', 'Count'])
cc5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_emcp.csv", sep=';', usecols=['SMILES', 'Count'])
cc6 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_emv.csv", sep=';', usecols=['SMILES', 'Count'])
cc7 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_sos.csv", sep=';', usecols=['SMILES', 'Count'])

gg1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_actcoeff.csv", sep=';', index_col= 'Unnamed: 0')
gg2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_dens.csv", sep=';', index_col= 'Unnamed: 0')
gg3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_dielec.csv", sep=';', index_col= 'Unnamed: 0')
gg4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_eme.csv", sep=';', index_col= 'Unnamed: 0')
gg5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_emcp.csv", sep=';', index_col= 'Unnamed: 0')
gg6 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_emv.csv", sep=';', index_col= 'Unnamed: 0')
gg7 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_sos.csv", sep=';', index_col= 'Unnamed: 0')

a = c2.merge(c1,how='outer',on=['SMILES'],suffixes=(' dens pure', ' cpmol pure')).merge(c5,how='outer',on=['SMILES'],suffixes=(' cpmol pure', ' sos pure')).merge(c3,how='outer',on=['SMILES'],suffixes=(' sos pure', ' dielec pure')).merge(c4,how='outer',on=['SMILES'],suffixes=(' dielec pure', ' hmol pure'))
a.replace(np.nan,0,inplace=True)
a.rename(columns={'Count':'Count hmol pure'},inplace=True)

b = cc2.merge(cc1,how='outer',on=['SMILES'],suffixes=(' dens binary', ' actcoeff binary')).merge(cc7,how='outer',on=['SMILES'],suffixes=(' actcoeff binary', ' sos binary')).merge(cc3,how='outer',on=['SMILES'],suffixes=(' sos binary', ' dielec binary')).merge(cc4,how='outer',on=['SMILES'],suffixes=(' dielec binary', ' eme binary')).merge(cc5,how='outer',on=['SMILES'],suffixes=(' eme binary', ' emcp binary')).merge(cc6,how='outer',on=['SMILES'],suffixes=(' emcp binary', ' emv binary'))
b.replace(np.nan,0,inplace=True)
b.rename(columns={'Count':'Count emv binary'},inplace=True)

bb = ccc2.merge(ccc1,how='outer',on=['SMILES'],suffixes=(' dens binary', ' actcoeff binary')).merge(ccc7,how='outer',on=['SMILES'],suffixes=(' actcoeff binary', ' sos binary')).merge(ccc3,how='outer',on=['SMILES'],suffixes=(' sos binary', ' dielec binary')).merge(ccc4,how='outer',on=['SMILES'],suffixes=(' dielec binary', ' eme binary')).merge(ccc5,how='outer',on=['SMILES'],suffixes=(' eme binary', ' emcp binary')).merge(ccc6,how='outer',on=['SMILES'],suffixes=(' emcp binary', ' emv binary'))
bb.replace(np.nan,0,inplace=True)
bb.rename(columns={'Count':'Count emv binary'},inplace=True)

c = pd.merge(a,b,how='outer',on=['SMILES'])
c.replace(np.nan,0,inplace=True)

cc = pd.merge(a,bb,how='outer',on=['SMILES'])
cc.replace(np.nan,0,inplace=True)
####################################################################################################################
a.insert(0,'Component',a.SMILES.apply(lambda x: resolve_cached(x, "iupac_name")))
cols = list(a)
cols.insert(1, cols.pop(cols.index('SMILES')))
a = a.ix[:, cols]
b.insert(0,'Component',b.SMILES.apply(lambda x: resolve_cached(x, "iupac_name")))
cols = list(b)
cols.insert(1, cols.pop(cols.index('SMILES')))
b = b.ix[:, cols]
bb.insert(0,'Component',bb.SMILES.apply(lambda x: resolve_cached(x, "iupac_name")))
cols = list(bb)
cols.insert(1, cols.pop(cols.index('SMILES')))
bb = bb.ix[:, cols]
c.insert(0,'Component',c.SMILES.apply(lambda x: resolve_cached(x, "iupac_name")))
cols = list(c)
cols.insert(1, cols.pop(cols.index('SMILES')))
c = c.ix[:, cols]
cc.insert(0,'Component',cc.SMILES.apply(lambda x: resolve_cached(x, "iupac_name")))
cols = list(cc)
cols.insert(1, cols.pop(cols.index('SMILES')))
cc = cc.ix[:, cols]

d = gg2.merge(gg1,how='outer',on=['Mixture'],suffixes=(' dens binary', ' actcoeff binary')).merge(gg7,how='outer',on=['Mixture'],suffixes=(' actcoeff binary', ' sos binary')).merge(gg3,how='outer',on=['Mixture'],suffixes=(' sos binary', ' dielec binary')).merge(gg4,how='outer',on=['Mixture'],suffixes=(' dielec binary', ' eme binary')).merge(gg5,how='outer',on=['Mixture'],suffixes=(' eme binary', ' emcp binary')).merge(gg6,how='outer',on=['Mixture'],suffixes=(' emcp binary', ' emv binary'))
d.replace(np.nan,0,inplace=True)
d.rename(columns={'Count':'Mixture Count emv binary'},inplace=True)


d["x1"], d["x2"] =  zip(*d["Mixture"].str.split('__').tolist())

cols = list(d)
cols.insert(0, cols.pop(cols.index('Mixture')))
cols.insert(1, cols.pop(cols.index('x1')))
cols.insert(2, cols.pop(cols.index('x2')))
d = d.ix[:, cols]

d.insert(3,'SMILES1',d.x1.apply(lambda x: resolve_cached(x, "smiles")))
d.insert(4,'SMILES2',d.x2.apply(lambda x: resolve_cached(x, "smiles")))

ind = np.logical_not(d.SMILES1.isin(S.SMILES) & d.SMILES2.isin(S.SMILES))

d_int = d[ind]

# Create total count column for each df by summing across all properties
a['Count all pure'] = a.sum(axis=1)
b['Count all binary'] = b.sum(axis=1)
bb['Count all binary'] = bb.sum(axis=1)
c['Count all'] = c.sum(axis=1)
cc['Count all'] = cc.sum(axis=1)
d['Mixture Count all'] = d.sum(axis=1)
d_int['Mixture Count all'] = d_int.sum(axis=1)

# Sort data frame by total count column in descending order 
a = a.sort_values(by='Count all pure', ascending=False)
b = b.sort_values(by='Count all binary', ascending=False)
bb = bb.sort_values(by='Count all binary', ascending=False)
c = c.sort_values(by='Count all', ascending=False)
cc = cc.sort_values(by='Count all', ascending=False)
d = d.sort_values(by='Mixture Count all', ascending=False)
d_int = d_int.sort_values(by='Mixture Count all', ascending=False)

# Save df to csv and pickle file format
a.to_csv("/home/bmanubay/.thermoml/tables/Ken/purecomp_counts_all.csv", sep=';')
b.to_csv("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_all.csv", sep=';')
bb.to_csv("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_interesting.csv", sep=';')
c.to_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_all.csv", sep=';')
cc.to_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_interesting.csv", sep=';')
d.to_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_all.csv", sep=';')
d_int.to_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_interesting.csv", sep=';')

a.to_pickle("/home/bmanubay/.thermoml/tables/Ken/purecomp_counts_all.pkl")
b.to_pickle("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_all.pkl")
bb.to_pickle("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_interesting.pkl")
c.to_pickle("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_all.pkl")
cc.to_pickle("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_interesting.pkl")
d.to_pickle("/home/bmanubay/.thermoml/tables/Ken/mix_counts_all.pkl")
d_int.to_pickle("/home/bmanubay/.thermoml/tables/Ken/mix_counts_interesting.pkl")