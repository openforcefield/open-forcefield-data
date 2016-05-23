# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:35:57 2016

@author: bmanubay
"""

import thermopyl as th 
from thermopyl import thermoml_lib
import cirpy
import numpy as np
import pandas as pd
from sklearn.externals.joblib import Memory

mem = Memory(cachedir="/home/bmanubay/.thermoml/")

@mem.cache
def resolve_cached(x, rtype):
   return cirpy.resolve(x, rtype)
   
# Compounds of most interest as decide by David, Chris and Bryce   
davmollist = ['2,2,4-trimethylpentane', 'cycloheptane', 'diisopropylether', 'isopropyl ether', 'dimethoxymethane', '2,3-dimethylbutane', '2,2-dimethylbutane', '3-methylpentane', 'neohexane', '4-methyl-2-pentanol', '2-methyl-2-pentanol', '1,1-diethoxyethane', 'tert-butanol', 'tetrahydrofuran', 'heptane', 'water', 'ethanol', '1-butanol', 'methyl tert-butyl ether']
S = pd.DataFrame({'IUPAC_Names': davmollist}, columns = ['IUPAC_Names'])
S["SMILES"] = S.IUPAC_Names.apply(lambda x: resolve_cached(x, "smiles"))

df = th.pandas_dataframe()
dt = list(df.columns)

bad_filenames = ["/home/bmanubay/.thermoml/j.fluid.2013.12.014.xml"]  # This file confirmed to have possible data entry errors.
df = df[~df.filename.isin(bad_filenames)]

experiments = ["Mass density, kg/m3","Speed of sound, m/s", "Relative permittivity at zero frequency", "Molar heat capacity at constant pressure, J/K/mol", "Molar enthalpy, kJ/mol"]

ind_list = [df[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
df = df.ix[ind]

name_to_formula = pd.read_hdf("/home/bmanubay/.thermoml/compound_name_to_formula.h5", 'data')
name_to_formula = name_to_formula.dropna()

# Extract rows with two components
df["n_components"] = df.components.apply(lambda x: len(x.split("__")))
df = df[df.n_components == 1]
df.dropna(axis=1, how='all', inplace=True)


df["formula"] = df.components.apply(lambda chemical: name_to_formula[chemical])

heavy_atoms = ["C", "O"]
desired_atoms = ["H"] + heavy_atoms

df["n_atoms"] = df.formula.apply(lambda formula_string : thermoml_lib.count_atoms(formula_string))
df["n_heavy_atoms"] = df.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, heavy_atoms))
df["n_desired_atoms"] = df.formula.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, desired_atoms))
df["n_other_atoms"] = df.n_atoms - df.n_desired_atoms

df = df[df.n_other_atoms == 0]

df = df[df.n_heavy_atoms > 0]
df = df[df.n_heavy_atoms <= 10]
df.dropna(axis=1, how='all', inplace=True)

df["SMILES"] = df.components.apply(lambda x: resolve_cached(x, "smiles"))  # This should be cached via sklearn.
df = df[df.SMILES != None]
df = df[df["SMILES"].str.contains('=O') == False] # Getting rid of data sets with C=O and C=C occurrences
df = df[df["SMILES"].str.contains('#') == False]
df = df[df["SMILES"].str.contains('O=') == False]
df = df[df["SMILES"].str.contains('=C') == False]
df = df[df["SMILES"].str.contains('C=') == False]
df.dropna(subset=["SMILES"], inplace=True)
df = df.ix[df.SMILES.dropna().index]

df["cas"] = df.components.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "cas")))  # This should be cached via sklearn.
df["InChI"] = df.components.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "stdinchikey")))
df = df[df.cas != None]
df = df.ix[df.cas.dropna().index]

# Neither names (components) nor smiles are unique.  Use CAS to ensure consistency.
cannonical_smiles_lookup = df.groupby("cas").SMILES.first()
cannonical_components_lookup = df.groupby("cas").components.first()


df["SMILES"] = df.cas.apply(lambda x: cannonical_smiles_lookup[x])
df["components"] = df.cas.apply(lambda x: cannonical_components_lookup[x])

# Extract rows with temperature between 128 and 399 K
df = df[df['Temperature, K'] > 250.]
df = df[df['Temperature, K'] < 400.]

# Extract rows with pressure between 101.325 kPa and 101325 kPa
df = df[df['Pressure, kPa'] > 100.]
df = df[df['Pressure, kPa'] < 102000.]

# Strip rows not in liquid phase
df = df[df['phase']=='Liquid']

df.dropna(axis=1, how='all', inplace=True)

df["filename"] = df["filename"].map(lambda x: x.lstrip('/home/bmanubay/.thermoml/').rstrip('.xml'))

def dfpretty(df, prop):
    dfbig = pd.concat([df['filename'], df["components"], df["SMILES"], df["cas"], df["InChI"], df["Temperature, K"], df["Pressure, kPa"], df[prop], df[prop+"_std"]], axis=1, keys=["filename", "components", "SMILES", "CAS", "InChI", "Temperature, K", "Pressure, kPa", prop, prop+"_std"])
    dfbig[prop+"_std"].replace('nan', np.nan, inplace=True)
    dfbig = dfbig[np.isnan(dfbig[prop+"_std"])==False]
    cannonical_smiles_lookup = dfbig.groupby("CAS").SMILES.first()
    cannonical_components_lookup = dfbig.groupby("CAS").components.first()
    dfbig["SMILES"] = dfbig.CAS.apply(lambda x: cannonical_smiles_lookup[x])
    dfbig["components"] = dfbig.CAS.apply(lambda x: cannonical_components_lookup[x])
    a = dfbig["filename"].value_counts()
    a = a.reset_index()
    a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
    b = dfbig["InChI"].value_counts()
    b = b.reset_index()
    b.rename(columns={"index":"InChI","InChI":"Count"},inplace=True)
    b["Component"] = b.InChI.apply(lambda x: resolve_cached(x, "iupac_name")) 
    b["SMILES"] = b.InChI.apply(lambda x: resolve_cached(x, "smiles")) 
    
    return dfbig, a, b
    
                   
dfbig = pd.concat([df['filename'], df["components"], df["SMILES"], df["cas"], df["InChI"], df["Temperature, K"], df["Pressure, kPa"], df["Mass density, kg/m3"], df["Mass density, kg/m3_std"], df["Speed of sound, m/s"], df["Speed of sound, m/s_std"], df["Relative permittivity at zero frequency"], df["Relative permittivity at zero frequency_std"], df["Molar heat capacity at constant pressure, J/K/mol"], df["Molar heat capacity at constant pressure, J/K/mol_std"], df["Molar enthalpy, kJ/mol"], df["Molar enthalpy, kJ/mol"]] , axis=1, keys=["filename", "components", "SMILES", "CAS", "InChI", "Temperature, K", "Pressure, kPa", "Mass density, kg/m3", "Mass density, kg/m3_std", "Speed of sound, m/s", "Speed of sound, m/s_std", "Relative permittivity at zero frequency", "Relative permittivity at zero frequency_std", "Molar heat capacity at constant pressure, J/K/mol", "Molar heat capacity at constant pressure, J/K/mol_std", "Molar enthalpy, kJ/mol", "Molar enthalpy, kJ/mol"])
a = dfbig["filename"].value_counts()
a = a.reset_index()
a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b = dfbig["InChI"].value_counts()    
b = b.reset_index()    
b.rename(columns={"index":"InChI","InChI":"Count"},inplace=True)    
b["Component"] = b.InChI.apply(lambda x: resolve_cached(x, "iupac_name"))     
b["SMILES"] = b.InChI.apply(lambda x: resolve_cached(x, "smiles"))

df1, a1, b1 = dfpretty(df, "Mass density, kg/m3")
df2, a2, b2 = dfpretty(df, "Speed of sound, m/s")
df3, a3, b3 = dfpretty(df, "Relative permittivity at zero frequency")
df4, a4, b4 = dfpretty(df, "Molar heat capacity at constant pressure, J/K/mol")
df5, a5, b5 = dfpretty(df, "Molar enthalpy, kJ/mol")


pathdf = "/home/bmanubay/.thermoml/tables/Ken/Pure/Property data/"
pathjourn = "/home/bmanubay/.thermoml/tables/Ken/Pure/Journal name counts/"
pathcomp = "/home/bmanubay/.thermoml/tables/Ken/Pure/Component counts/"

def saveprettycsv(df, path, filename):
    df.to_csv(path+filename, sep =';')

def saveprettypickle(df, path, filename):
    df.to_pickle(path+filename)

# save csv with ; delimiter
saveprettycsv(dfbig, pathdf, "alldata_pure.csv")    
saveprettycsv(df1, pathdf, "dens_pure.csv") 
saveprettycsv(df2, pathdf, "sos_pure.csv") 
saveprettycsv(df3, pathdf, "dielec_pure.csv") 
saveprettycsv(df4, pathdf, "cpmol_pure.csv") 
saveprettycsv(df5, pathdf, "hmol_pure.csv") 

saveprettycsv(a, pathjourn, "purename_counts_all.csv")
saveprettycsv(a1, pathjourn, "purename_counts_dens.csv")
saveprettycsv(a2, pathjourn, "purename_counts_sos.csv")
saveprettycsv(a3, pathjourn, "purename_counts_dielec.csv")
saveprettycsv(a4, pathjourn, "purename_counts_cpmol.csv")
saveprettycsv(a5, pathjourn, "purename_counts_hmol.csv")

saveprettycsv(b, pathcomp, "purecomp_counts_all.csv")
saveprettycsv(b1, pathcomp, "purecomp_counts_dens.csv")
saveprettycsv(b2, pathcomp, "purecomp_counts_sos.csv")
saveprettycsv(b3, pathcomp, "purecomp_counts_dielec.csv")
saveprettycsv(b4, pathcomp, "purecomp_counts_cpmol.csv")
saveprettycsv(b5, pathcomp, "purecomp_counts_hmol.csv")

# save pickle
saveprettypickle(dfbig, pathdf, "alldata_pure.pkl")    
saveprettypickle(df1, pathdf, "dens_pure.pkl") 
saveprettypickle(df2, pathdf, "sos_pure.pkl") 
saveprettypickle(df3, pathdf, "dielec_pure.pkl") 
saveprettypickle(df4, pathdf, "cpmol_pure.pkl") 
saveprettypickle(df5, pathdf, "hmol_pure.pkl") 

saveprettypickle(a, pathjourn, "purename_counts_all.pkl")
saveprettypickle(a1, pathjourn, "purename_counts_dens.pkl")
saveprettypickle(a2, pathjourn, "purename_counts_sos.pkl")
saveprettypickle(a3, pathjourn, "purename_counts_dielec.pkl")
saveprettypickle(a4, pathjourn, "purename_counts_cpmol.pkl")
saveprettypickle(a5, pathjourn, "purename_counts_hmol.pkl")

saveprettypickle(b, pathcomp, "purecomp_counts_all.pkl")
saveprettypickle(b1, pathcomp, "purecomp_counts_dens.pkl")
saveprettypickle(b2, pathcomp, "purecomp_counts_sos.pkl")
saveprettypickle(b3, pathcomp, "purecomp_counts_dielec.pkl")
saveprettypickle(b4, pathcomp, "purecomp_counts_cpmol.pkl")
saveprettypickle(b5, pathcomp, "purecomp_counts_hmol.pkl") 


def SMILESchk(df, SMILES):
    ind = df.SMILES.isin(SMILES)
    dfDavSet = df[ind]
    
    return dfDavSet

dfD = SMILESchk(df, S.SMILES)
df1D, a1D, b1D = dfpretty(dfD, "Mass density, kg/m3")
df2D, a2D, b2D = dfpretty(dfD, "Speed of sound, m/s")
df3D, a3D, b3D = dfpretty(dfD, "Relative permittivity at zero frequency")
df4D, a4D, b4D = dfpretty(dfD, "Molar heat capacity at constant pressure, J/K/mol")
df5D, a5D, b5D = dfpretty(dfD, "Molar enthalpy, kJ/mol")

# save csv with ; delimiter  
saveprettycsv(df1D, pathdf, "dens_pure_BDC.csv") 
saveprettycsv(df2D, pathdf, "sos_pure_BDC.csv") 
saveprettycsv(df3D, pathdf, "dielec_pure_BDC.csv") 
saveprettycsv(df4D, pathdf, "cpmol_pure_BDC.csv") 
saveprettycsv(df5D, pathdf, "hmol_pure_BDC.csv") 

saveprettycsv(a1D, pathjourn, "purename_counts_dens_BDC.csv")
saveprettycsv(a2D, pathjourn, "purename_counts_sos_BDC.csv")
saveprettycsv(a3D, pathjourn, "purename_counts_dielec_BDC.csv")
saveprettycsv(a4D, pathjourn, "purename_counts_cpmol_BDC.csv")
saveprettycsv(a5D, pathjourn, "purename_counts_hmol_BDC.csv")

saveprettycsv(b1D, pathcomp, "purecomp_counts_dens_BDC.csv")
saveprettycsv(b2D, pathcomp, "purecomp_counts_sos_BDC.csv")
saveprettycsv(b3D, pathcomp, "purecomp_counts_dielec_BDC.csv")
saveprettycsv(b4D, pathcomp, "purecomp_counts_cpmol_BDC.csv")
saveprettycsv(b5D, pathcomp, "purecomp_counts_hmol_BDC.csv")

# save pickle 
saveprettypickle(df1D, pathdf, "dens_pure_BDC.pkl") 
saveprettypickle(df2D, pathdf, "sos_pure_BDC.pkl") 
saveprettypickle(df3D, pathdf, "dielec_pure_BDC.pkl") 
saveprettypickle(df4D, pathdf, "cpmol_pure_BDC.pkl") 
saveprettypickle(df5D, pathdf, "hmol_pure_BDC.pkl") 

saveprettypickle(a1D, pathjourn, "purename_counts_dens_BDC.pkl")
saveprettypickle(a2D, pathjourn, "purename_counts_sos_BDC.pkl")
saveprettypickle(a3D, pathjourn, "purename_counts_dielec_BDC.pkl")
saveprettypickle(a4D, pathjourn, "purename_counts_cpmol_BDC.pkl")
saveprettypickle(a5D, pathjourn, "purename_counts_hmol_BDC.pkl")

saveprettypickle(b1D, pathcomp, "purecomp_counts_dens_BDC.pkl")
saveprettypickle(b2D, pathcomp, "purecomp_counts_sos_BDC.pkl")
saveprettypickle(b3D, pathcomp, "purecomp_counts_dielec_BDC.pkl")
saveprettypickle(b4D, pathcomp, "purecomp_counts_cpmol_BDC.pkl")
saveprettypickle(b5D, pathcomp, "purecomp_counts_hmol_BDC.pkl") 


