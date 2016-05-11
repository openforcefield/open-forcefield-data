# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:52:43 2016

@author: bmanubay
"""

# Check what mixtures we have data on after eliminating alkane-alkane mixes
import cirpy
import numpy as np
import pandas as pd
from sklearn.externals.joblib import Memory

mem = Memory(cachedir="/home/bmanubay/.thermoml/")

@mem.cache
def resolve_cached(x, rtype):
   return cirpy.resolve(x, rtype)

# Read in mix_counts_all.csv
a = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_all.csv")
b = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_chain_filt1.smi",delim_whitespace=True,header=None)
c = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_rings_filt1.smi",delim_whitespace=True,header=None)

d = pd.merge(b,c,on=0,how='outer')
d.columns = ["SMILES","chains,","rings"]

# Define list of all alkane SMILES strings that appear in allcomp_counts_all.csv
SMILES_alk = ['C', 'CC', 'CCC', 'CCCC', 'CCCCC', 'CCCCCC', 'CCCCCCC', 'CCCCCCCC', 'CCCCCCCCC', 'CCCCCCCCCC', 'CC(C)C', 'CCC(C)C', 'CCCC(C)C', 'C1CCCCC1', 'CC1CCCCC1', 'CCCCCC(C)C', 'CC(C)C(C)C', 'CCC(C)(C)C', 'CCC(C)CC', 'CCCC(C)C', 'CC(C)CC(C)(C)C', 'C1CCCC1', 'C1CCCCCCC1', 'CCC1CCCCC1', 'CC1CCC(C)C(C)C1', 'CCCCC1CCCCC1', 'CC1CCCC1', 'CCCCCCC(C)C', 'CCCCCCCC(C)C', 'CCCCC(C)CCC', 'CCC(C)CCC(C)CC', 'CCC(C)CC(C)CC', 'CCCCCC(C)C', 'C1CCCCCC1', 'CC(C)C(C)C', 'CCC(C)(C)C']

S = pd.DataFrame({'SMILES': SMILES_alk}, columns = ['SMILES'])

a["x1"], a["x2"] =  zip(*a["Mixture"].str.split('__').tolist())

a["SMILES1"] = a.x1.apply(lambda x: resolve_cached(x, "smiles"))
a["SMILES2"] = a.x2.apply(lambda x: resolve_cached(x, "smiles"))

ind = a.SMILES1.isin(d.SMILES) & a.SMILES2.isin(d.SMILES)

a = a[ind]

ind1 = np.logical_not(a.SMILES1.isin(S.SMILES) & a.SMILES2.isin(S.SMILES))

a_mix_chk = a[ind1]

a_mix_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_interesting.csv") 
