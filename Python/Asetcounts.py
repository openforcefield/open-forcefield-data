# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:16:31 2016

@author: bmanubay
"""

import thermopyl as th 
from thermopyl import thermoml_lib
import cirpy
import numpy as np
import pandas as pd

mem = Memory(cachedir="/home/bmanubay/.thermoml/")

@mem.cache
def resolve_cached(x, rtype):
   return cirpy.resolve(x, rtype)


# Load in hand filtered property csv corresponding to A set

## Pure compounds 
a1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/cpmol_pure_A_set.csv")
a2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/dens_pure_A_set.csv")
a3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/dielec_pure_A_set.csv")
a4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/hmol_pure_A_set.csv")
a5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/sos_pure_A_set.csv")
a6 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/vmol_pure_A_set.csv")
a7 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/vspec_pure_A_set.csv")

## Binary mixtures
aa1 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/actcoeff_bin_A_set.csv")
aa2 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/dens_bin_A_set.csv")
aa3 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/dielec_bin_A_set.csv")
aa4 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/eme_bin_A_set.csv")
aa5 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/emhp_bin_A_set.csv")
aa6 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/emv_bin_A_set.csv")
aa7 = pd.read_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/sos_bin_A_set.csv")

# Get counts based on journal name

## Pure compounds 
b1 = a1["filename"].value_counts()
b1 = b1.reset_index()
b1.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b2 = a2["filename"].value_counts()
b2 = b2.reset_index()
b2.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b3 = a3["filename"].value_counts()
b3 = b3.reset_index()
b3.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b4 = a4["filename"].value_counts()
b4 = b4.reset_index()
b4.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b5 = a5["filename"].value_counts()
b5 = b5.reset_index()
b5.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b6 = a6["filename"].value_counts()
b6 = b6.reset_index()
b6.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b7 = a7["filename"].value_counts()
b7 = b7.reset_index()
b7.rename(columns={"index":"Filename","filename":"Count"},inplace=True)

## Binary mixtures
bb1 = aa1["filename"].value_counts()
bb1 = bb1.reset_index()
bb1.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bb2 = aa2["filename"].value_counts()
bb2 = bb2.reset_index()
bb2.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bb3 = aa3["filename"].value_counts()
bb3 = bb3.reset_index()
bb3.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bb4 = aa4["filename"].value_counts()
bb4 = bb4.reset_index()
bb4.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bb5 = aa5["filename"].value_counts()
bb5 = bb5.reset_index()
bb5.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bb6 = aa6["filename"].value_counts()
bb6 = bb6.reset_index()
bb6.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bb7 = aa7["filename"].value_counts()
bb7 = bb7.reset_index()
bb7.rename(columns={"index":"Filename","filename":"Count"},inplace=True)

# Get counts based on component names

## Pure compounds
c1 = a1["components"].value_counts()
c1 = c1.reset_index()
c1.rename(columns={"index":"Component","components":"Count"},inplace=True)
c1 = pd.concat([c1["Component"], c1["Count"]], axis=1, keys=["Component", "Count"]) 
c2 = a2["components"].value_counts()
c2 = c2.reset_index()
c2.rename(columns={"index":"Component","components":"Count"},inplace=True)
c2 = pd.concat([c2["Component"], c2["Count"]], axis=1, keys=["Component", "Count"]) 
c3 = a3["components"].value_counts()
c3 = c3.reset_index()
c3.rename(columns={"index":"Component","components":"Count"},inplace=True)
c3 = pd.concat([c3["Component"], c3["Count"]], axis=1, keys=["Component", "Count"])
c4 = a4["components"].value_counts()
c4 = c4.reset_index()
c4.rename(columns={"index":"Component","components":"Count"},inplace=True)
c4 = pd.concat([c4["Component"], c4["Count"]], axis=1, keys=["Component", "Count"])  
c5 = a5["components"].value_counts()
c5 = c5.reset_index()
c5.rename(columns={"index":"Component","components":"Count"},inplace=True)
c5 = pd.concat([c5["Component"], c5["Count"]], axis=1, keys=["Component", "Count"]) 
c6 = a6["components"].value_counts()
c6 = c6.reset_index()
c6.rename(columns={"index":"Component","components":"Count"},inplace=True)
c6 = pd.concat([c6["Component"], c6["Count"]], axis=1, keys=["Component", "Count"]) 
c7 = a7["components"].value_counts()
c7 = c7.reset_index()
c7.rename(columns={"index":"Component","components":"Count"},inplace=True)
c7 = pd.concat([c7["Component"], c7["Count"]], axis=1, keys=["Component", "Count"]) 

## Binary Mixtures
cc1c = pd.concat([aa1["x1"], aa1["x2"]], axis=1, keys = ["x1", "x2"]) 
count1c = pd.Series(cc1c.squeeze().values.ravel()).value_counts()
count1c = count1c.reset_index()
count1c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc1 = pd.concat([count1c["Component"], count1c["Count"]], axis=1, keys=["Component", "Component Count"])
cc2c = pd.concat([aa2["x1"], aa2["x2"]], axis=1, keys = ["x1", "x2"]) 
count2c = pd.Series(cc2c.squeeze().values.ravel()).value_counts()
count2c = count2c.reset_index()
count2c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc2 = pd.concat([count2c["Component"], count2c["Count"]], axis=1, keys=["Component", "Component Count"])
cc3c = pd.concat([aa3["x1"], aa3["x2"]], axis=1, keys = ["x1", "x2"]) 
count3c = pd.Series(cc3c.squeeze().values.ravel()).value_counts()
count3c = count3c.reset_index()
count3c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc3 = pd.concat([count3c["Component"], count3c["Count"]], axis=1, keys=["Component", "Component Count"])
cc4c = pd.concat([aa4["x1"], aa4["x2"]], axis=1, keys = ["x1", "x2"]) 
count4c = pd.Series(cc4c.squeeze().values.ravel()).value_counts()
count4c = count4c.reset_index()
count4c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc4 = pd.concat([count4c["Component"], count4c["Count"]], axis=1, keys=["Component", "Component Count"])
cc5c = pd.concat([aa5["x1"], aa5["x2"]], axis=1, keys = ["x1", "x2"]) 
count5c = pd.Series(cc5c.squeeze().values.ravel()).value_counts()
count5c = count5c.reset_index()
count5c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc5 = pd.concat([count5c["Component"], count5c["Count"]], axis=1, keys=["Component", "Component Count"])
cc6c = pd.concat([aa6["x1"], aa6["x2"]], axis=1, keys = ["x1", "x2"]) 
count6c = pd.Series(cc6c.squeeze().values.ravel()).value_counts()
count6c = count6c.reset_index()
count6c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc6 = pd.concat([count6c["Component"], count6c["Count"]], axis=1, keys=["Component", "Component Count"])
cc7c = pd.concat([aa7["x1"], aa7["x2"]], axis=1, keys = ["x1", "x2"]) 
count7c = pd.Series(cc7c.squeeze().values.ravel()).value_counts()
count7c = count7c.reset_index()
count7c.rename(columns={"index":"Component",0:"Count"},inplace=True)
cc7 = pd.concat([count7c["Component"], count7c["Count"]], axis=1, keys=["Component", "Component Count"])

# Get binary mixture counts
gg1 = aa1["components"].value_counts()
gg1 = gg1.reset_index()
gg1.rename(columns={"index":"Mixture","components":"Count"},inplace=True)
gg2 = aa2["components"].value_counts()
gg2 = gg2.reset_index()
gg2.rename(columns={"index":"Mixture","components":"Count"},inplace=True)
gg3 = aa3["components"].value_counts()
gg3 = gg3.reset_index()
gg3.rename(columns={"index":"Mixture","components":"Count"},inplace=True)
gg4 = aa4["components"].value_counts()
gg4 = gg4.reset_index()
gg4.rename(columns={"index":"Mixture","components":"Count"},inplace=True)
gg5 = aa5["components"].value_counts()
gg5 = gg5.reset_index()
gg5.rename(columns={"index":"Mixture","components":"Count"},inplace=True)
gg6 = aa6["components"].value_counts()
gg6 = gg6.reset_index()
gg6.rename(columns={"index":"Mixture","components":"Count"},inplace=True)
gg7 = aa7["components"].value_counts()
gg7 = gg7.reset_index()
gg7.rename(columns={"index":"Mixture","components":"Count"},inplace=True)


# Save pure journal counts as .csv
b1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_cpmol_A_set.csv")
b2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_dens_A_set.csv")
b3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_dielec_A_set.csv")
b4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_hmol_A_set.csv")
b5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_sos_A_set.csv")
b6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_vmol_A_set.csv")
b7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/purename_counts_vspec_A_set.csv")

# Save binary journal counts as .csv
bb1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_actcoeff_A_set.csv")
bb2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_dens_A_set.csv")
bb3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_dielec_A_set.csv")
bb4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_eme_A_set.csv")
bb5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_emhp_A_set.csv")
bb6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_emv_A_set.csv")
bb7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_sos_A_set.csv")

# Save pure component counts as .csv
c1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_cpmol_A_set.csv")
c2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_dens_A_set.csv")
c3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_dielec_A_set.csv")
c4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_hmol_A_set.csv")
c5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_sos_A_set.csv")
c6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_vmol_A_set.csv")
c7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/purecomp_counts_vspec_A_set.csv")

# Save binary component counts as .csv
cc1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_actcoeff_A_set.csv")
cc2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_dens_A_set.csv")
cc3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_dielec_A_set.csv")
cc4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_eme_A_set.csv")
cc5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_emhp_A_set.csv")
cc6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_emv_A_set.csv")
cc7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/bincomp_counts_sos_A_set.csv")

# Save binary mixture counts as .csv
gg1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_actcoeff_A_set.csv")
gg2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_dens_A_set.csv")
gg3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_dielec_A_set.csv")
gg4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_eme_A_set.csv")
gg5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_emhp_A_set.csv")
gg6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_emv_A_set.csv")
gg7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_sos_A_set.csv")









