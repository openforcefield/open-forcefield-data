# -*- coding: utf-8 -*-
"""
Created on Sun May  8 18:29:53 2016

@author: bmanubay
"""

# Check what moelcules we have appear in Chris's list

import pandas as pd

a = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_all.csv")
b = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_chain_filt1.smi",delim_whitespace=True,header=None)
c = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/AlkEthOH_distrib2016may01/AlkEthOH_rings_filt1.smi",delim_whitespace=True,header=None)

d = pd.merge(b,c,on=0,how='outer')
d.columns = ["SMILES","chains,","rings"]

ind = a.SMILES.isin(d.SMILES)

a_diverse_chk = a[ind]


a_diverse_chk.to_csv("/home/bmanubay/.thermoml/tables/Ken/allcomp_counts_diverse.csv") 