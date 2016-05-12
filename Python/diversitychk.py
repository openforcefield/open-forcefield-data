# -*- coding: utf-8 -*-
"""
Created on Sun May  8 18:29:53 2016

@author: bmanubay
"""

# Check what moelcules we have appear in Chris's list

import pandas as pd

a0 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_all.csv")
a1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_interesting.csv")
a2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_all.csv")
a3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_interesting.csv")
a4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_all.csv")
a5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_interesting.csv")
a6 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/purecomp_counts_all.csv")
b = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_chain_filt1.smi",delim_whitespace=True,header=None)
c = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_rings_filt1.smi",delim_whitespace=True,header=None)

d = pd.merge(b,c,on=0,how='outer')
d.columns = ["SMILES","chains,","rings"]

ind0 = a0.SMILES.isin(d.SMILES)
ind1 = a1.SMILES.isin(d.SMILES)
ind2 = a2.SMILES.isin(d.SMILES)
ind3 = a3.SMILES.isin(d.SMILES)
ind4 = a4.SMILES1.isin(d.SMILES) & a4.SMILES2.isin(d.SMILES) 
ind5 = a5.SMILES1.isin(d.SMILES) & a5.SMILES2.isin(d.SMILES) 
ind6 = a6.SMILES.isin(d.SMILES)

a0_diverse_chk = a0[ind0]
a0_diverse_chk = a0_diverse_chk.drop("Unnamed: 0", axis = 1)
a1_diverse_chk = a1[ind1]
a1_diverse_chk = a1_diverse_chk.drop("Unnamed: 0", axis = 1)
a2_diverse_chk = a2[ind2]
a2_diverse_chk = a2_diverse_chk.drop("Unnamed: 0", axis = 1)
a3_diverse_chk = a3[ind3]
a3_diverse_chk = a3_diverse_chk.drop("Unnamed: 0", axis = 1)
a4_diverse_chk = a4[ind4]
a4_diverse_chk = a4_diverse_chk.drop("Unnamed: 0", axis = 1)
a5_diverse_chk = a5[ind5]
a5_diverse_chk = a5_diverse_chk.drop("Unnamed: 0", axis = 1)
a6_diverse_chk = a6[ind6]
a6_diverse_chk = a6_diverse_chk.drop("Unnamed: 0", axis = 1)

a0_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_diverse.csv") 
a1_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_interesting_diverse.csv")
a2_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_diverse.csv")
a3_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/bincomp_counts_interesting_diverse.csv")
a4_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_diverse.csv")
a5_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/mix_counts_intersting_diverse.csv")
a6_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/purecomp_counts_diverse.csv")