# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:33:13 2016

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

experiments = ["Mass density, kg/m3", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol", "Excess molar heat capacity, J/K/mol", "Excess molar volume, m3/mol", "Activity coefficient", "Speed of sound, m/s", "Relative permittivity at zero frequency"]

ind_list = [df[exp].dropna().index for exp in experiments]
ind = reduce(lambda x,y: x.union(y), ind_list)
df = df.ix[ind]

name_to_formula = pd.read_hdf("/home/bmanubay/.thermoml/compound_name_to_formula.h5", 'data')
name_to_formula = name_to_formula.dropna()

# Extract rows with two components
df["n_components"] = df.components.apply(lambda x: len(x.split("__")))
df = df[df.n_components == 2]
df.dropna(axis=1, how='all', inplace=True)

counts_data = {}
counts_data["0.  Two Components"] = df.count()[experiments]

# Split components into separate columns (to use name_to_formula)
df["x1"], df["x2"] =  zip(*df["components"].str.split('__').tolist())
df['x2'].replace('', np.nan, inplace=True)
df.dropna(subset=['x2'], inplace=True)

df["formula1"] = df.x1.apply(lambda chemical: name_to_formula[chemical])
df["formula2"] = df.x2.apply(lambda chemical: name_to_formula[chemical])

heavy_atoms = ["C", "O"]
desired_atoms = ["H"] + heavy_atoms

df["n_atoms1"] = df.formula1.apply(lambda formula_string : thermoml_lib.count_atoms(formula_string))
df["n_heavy_atoms1"] = df.formula1.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, heavy_atoms))
df["n_desired_atoms1"] = df.formula1.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, desired_atoms))
df["n_other_atoms1"] = df.n_atoms1 - df.n_desired_atoms1
df["n_atoms2"] = df.formula2.apply(lambda formula_string : thermoml_lib.count_atoms(formula_string))
df["n_heavy_atoms2"] = df.formula2.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, heavy_atoms))
df["n_desired_atoms2"] = df.formula2.apply(lambda formula_string : thermoml_lib.count_atoms_in_set(formula_string, desired_atoms))
df["n_other_atoms2"] = df.n_atoms2 - df.n_desired_atoms2

df = df[df.n_other_atoms1 == 0]
df = df[df.n_other_atoms2 == 0]

counts_data["1.  Druglike Elements"] = df.count()[experiments]

df = df[df.n_heavy_atoms1 > 0]
df = df[df.n_heavy_atoms1 <= 10]
df = df[df.n_heavy_atoms2 > 0]
df = df[df.n_heavy_atoms2 <= 10]
df.dropna(axis=1, how='all', inplace=True)

counts_data["2.  Heavy Atoms"] = df.count()[experiments]

df["SMILES1"] = df.x1.apply(lambda x: resolve_cached(x, "smiles"))  # This should be cached via sklearn.
df = df[df.SMILES1 != None]
df = df[df["SMILES1"].str.contains('=O') == False] # Getting rid of data sets with C=O and C=C occurrences
df = df[df["SMILES1"].str.contains('#') == False]
df = df[df["SMILES1"].str.contains('O=') == False]
df = df[df["SMILES1"].str.contains('=C') == False]
df = df[df["SMILES1"].str.contains('C=') == False]
df.dropna(subset=["SMILES1"], inplace=True)
df = df.ix[df.SMILES1.dropna().index]
df["SMILES2"] = df.x2.apply(lambda x: resolve_cached(x, "smiles"))  # This should be cached via sklearn.
df = df[df.SMILES2 != None]
df = df[df["SMILES2"].str.contains('=O') == False] # Getting rid of data sets with C=O and C=C occurrences
df = df[df["SMILES2"].str.contains('#') == False]
df = df[df["SMILES2"].str.contains('O=') == False]
df = df[df["SMILES2"].str.contains('=C') == False]
df = df[df["SMILES2"].str.contains('C=') == False]
df.dropna(subset=["SMILES2"], inplace=True)
df = df.ix[df.SMILES2.dropna().index]

    
df["cas1"] = df.x1.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "cas")))  # This should be cached via sklearn.
df["InChI1"] = df.x1.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "stdinchikey")))
df = df[df.cas1 != None]
df = df.ix[df.cas1.dropna().index]
df["cas2"] = df.x2.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "cas")))  # This should be cached via sklearn.
df["InChI2"] = df.x2.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "stdinchikey")))
df = df[df.cas2 != None]
df = df.ix[df.cas2.dropna().index]


# Neither names (components) nor smiles are unique.  Use CAS to ensure consistency.
cannonical_smiles_lookup1 = df.groupby("cas1").SMILES1.first()
cannonical_components_lookup1 = df.groupby("cas1").x1.first()
cannonical_smiles_lookup2 = df.groupby("cas2").SMILES2.first()
cannonical_components_lookup2 = df.groupby("cas2").x2.first()


df["SMILES1"] = df.cas1.apply(lambda x: cannonical_smiles_lookup1[x])
df["x1"] = df.cas1.apply(lambda x: cannonical_components_lookup1[x])
df["SMILES2"] = df.cas2.apply(lambda x: cannonical_smiles_lookup2[x])
df["x2"] = df.cas2.apply(lambda x: cannonical_components_lookup2[x])

# Extract rows with temperature between 128 and 399 K
df = df[df['Temperature, K'] > 250.]
df = df[df['Temperature, K'] < 400.]

counts_data["3.  Temperature"] = df.count()[experiments]

# Extract rows with pressure between 101.325 kPa and 101325 kPa
df = df[df['Pressure, kPa'] > 100.]
df = df[df['Pressure, kPa'] < 102000.]

counts_data["4.  Pressure"] = df.count()[experiments]

# Strip rows not in liquid phase
df = df[df['phase']=='Liquid']

counts_data["5.  Liquid state"] = df.count()[experiments]


df.dropna(axis=1, how='all', inplace=True)

df["filename"] = df["filename"].map(lambda x: x.lstrip('/home/bmanubay/.thermoml/').rstrip('.xml'))

def dfpretty(df, prop):
    dfbig = pd.concat([df['filename'], df["x1"], df['x2'], df["SMILES1"], df["SMILES2"], df["cas1"], df["cas2"], df["InChI1"], df["InChI2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df[prop], df[prop+"_std"]], axis=1, keys=["filename", "x1", "x2", "SMILES1", "SMILES2", "cas1", "cas2", "InChI1", "InChI2", "components", "Mole fraction", "Temperature, K", "Pressure, kPa", prop, prop+"_std"])
    dfbig[prop+"_std"].replace('nan', np.nan, inplace=True)
    dfbig = dfbig[np.isnan(dfbig[prop+"_std"])==False]
    cannonical_smiles_lookup1 = dfbig.groupby("cas1").SMILES1.first()
    cannonical_components_lookup1 = dfbig.groupby("cas1").x1.first()
    cannonical_smiles_lookup2 = dfbig.groupby("cas2").SMILES2.first()
    cannonical_components_lookup2 = dfbig.groupby("cas2").x2.first()
    dfbig["SMILES1"] = dfbig.cas1.apply(lambda x: cannonical_smiles_lookup1[x])
    dfbig["x1"] = dfbig.cas1.apply(lambda x: cannonical_components_lookup1[x])
    dfbig["SMILES2"] = dfbig.cas2.apply(lambda x: cannonical_smiles_lookup2[x])
    dfbig["x2"] = dfbig.cas2.apply(lambda x: cannonical_components_lookup2[x])
    a = dfbig["filename"].value_counts()
    a = a.reset_index()
    a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
    bInChI = pd.concat([dfbig["InChI1"], dfbig["InChI2"]], axis=1, keys = ["InChI1", "InChI2"])
    b = pd.Series(bInChI.squeeze().values.ravel()).value_counts()
    b = b.reset_index()
    b.rename(columns={"index":"InChI",0:"Count"},inplace=True)
    b["Component"] = b.InChI.apply(lambda x: resolve_cached(x, "iupac_name"))
    b["SMILES"] = b.InChI.apply(lambda x: resolve_cached(x, "smiles"))
    c = dfbig["components"].value_counts()
    c = c.reset_index()
    c.rename(columns={"index":"Mixture", "components":"Count"}, inplace=True)
    
    return dfbig, a, b, c
    
    
dfbig = pd.concat([df['filename'], df['x1'], df['x2'], df["SMILES1"], df["SMILES2"], df["cas1"], df["cas2"], df["InChI1"], df["InChI2"], df["components"], df["Mole fraction"], df["Temperature, K"], df["Pressure, kPa"], df["Mass density, kg/m3"], df["Mass density, kg/m3_std"], df["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol"], df["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std"], df["Excess molar heat capacity, J/K/mol"], df["Excess molar heat capacity, J/K/mol_std"], df["Excess molar volume, m3/mol"], df["Excess molar volume, m3/mol_std"], df["Activity coefficient"], df["Activity coefficient_std"], df["Speed of sound, m/s"], df["Speed of sound, m/s_std"], df["Relative permittivity at zero frequency"], df["Relative permittivity at zero frequency_std"]], axis=1, keys=["filename", "x1", "x2", "SMILES1", "SMILES2", "cas1", "cas2", "InChI1", "InChI2", "components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Mass density, kg/m3", "Mass density, kg/m3_std", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std", "Excess molar heat capacity, J/K/mol", "Excess molar heat capacity, J/K/mol_std", "Excess molar volume, m3/mol", "Excess molar volume, m3/mol_std", "Activity coefficient", "Activity coefficient_std", "Speed of sound, m/s", "Speed of sound, m/s_std", "Relative permittivity at zero frequency", "Relative permittivity at zero frequency_std"])
cannonical_smiles_lookup1 = dfbig.groupby("cas1").SMILES1.first()
cannonical_components_lookup1 = dfbig.groupby("cas1").x1.first()
cannonical_smiles_lookup2 = dfbig.groupby("cas2").SMILES2.first()
cannonical_components_lookup2 = dfbig.groupby("cas2").x2.first()
dfbig["SMILES1"] = dfbig.cas1.apply(lambda x: cannonical_smiles_lookup1[x])
dfbig["x1"] = dfbig.cas1.apply(lambda x: cannonical_components_lookup1[x])
dfbig["SMILES2"] = dfbig.cas2.apply(lambda x: cannonical_smiles_lookup2[x])
dfbig["x2"] = dfbig.cas2.apply(lambda x: cannonical_components_lookup2[x])
a = dfbig["filename"].value_counts()
a = a.reset_index()
a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
bInChI = pd.concat([dfbig["InChI1"], dfbig["InChI2"]], axis=1, keys = ["InChI1", "InChI2"])
b = pd.Series(bInChI.squeeze().values.ravel()).value_counts()
b = b.reset_index()
b.rename(columns={"index":"InChI",0:"Count"},inplace=True)
b["Component"] = b.InChI.apply(lambda x: resolve_cached(x, "iupac_name"))
b["SMILES"] = b.InChI.apply(lambda x: resolve_cached(x, "smiles"))
c = dfbig["components"].value_counts()
c = c.reset_index()
c.rename(columns={"index":"Mixture", "components":"Count"}, inplace=True)    

df1, a1, b1, c1 = dfpretty(df, "Mass density, kg/m3")
df2, a2, b2, c2 = dfpretty(df, "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol")
df3, a3, b3, c3 = dfpretty(df, "Excess molar heat capacity, J/K/mol")
df4, a4, b4, c4 = dfpretty(df, "Excess molar volume, m3/mol")
df5, a5, b5, c5 = dfpretty(df, "Activity coefficient")
df6, a6, b6, c6 = dfpretty(df, "Speed of sound, m/s")
df7, a7, b7, c7 = dfpretty(df, "Relative permittivity at zero frequency")


pathdf = "/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/"
pathjourn = "/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/"
pathcomp = "/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/"
pathmix = "/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/"

def saveprettycsv(df, path, filename):
    df.to_csv(path+filename, sep =';')

def saveprettypickle(df, path, filename):
    df.to_pickle(path+filename)

# save csv with ; delimiter
saveprettycsv(dfbig, pathdf, "alldata_bin.csv")
saveprettycsv(df1, pathdf, "dens_bin.csv")
saveprettycsv(df2, pathdf, "eme_bin.csv")
saveprettycsv(df3, pathdf, "emcp_bin.csv")
saveprettycsv(df4, pathdf, "emv_bin.csv")
saveprettycsv(df5, pathdf, "actcoeff_bin.csv")
saveprettycsv(df6, pathdf, "sos_bin.csv")
saveprettycsv(df7, pathdf, "dielec_bin.csv")

saveprettycsv(a, pathjourn, "binname_counts_all.csv")
saveprettycsv(a1, pathjourn, "binname_counts_dens.csv")
saveprettycsv(a2, pathjourn, "binname_counts_eme.csv")
saveprettycsv(a3, pathjourn, "binname_counts_emcp.csv")
saveprettycsv(a4, pathjourn, "binname_counts_emv.csv")
saveprettycsv(a5, pathjourn, "binname_counts_actcoeff.csv")
saveprettycsv(a6, pathjourn, "binname_counts_sos.csv")
saveprettycsv(a7, pathjourn, "binname_counts_dielec.csv")

saveprettycsv(b, pathcomp, "bincomp_counts_all.csv")
saveprettycsv(b1, pathcomp, "bincomp_counts_dens.csv")
saveprettycsv(b2, pathcomp, "bincomp_counts_eme.csv")
saveprettycsv(b3, pathcomp, "bincomp_counts_emcp.csv")
saveprettycsv(b4, pathcomp, "bincomp_counts_emv.csv")
saveprettycsv(b5, pathcomp, "bincomp_counts_actcoeff.csv")
saveprettycsv(b6, pathcomp, "bincomp_counts_sos.csv")
saveprettycsv(b7, pathcomp, "bincomp_counts_dielec.csv")

saveprettycsv(c, pathmix, "mix_counts_all.csv")
saveprettycsv(c1, pathmix, "mix_counts_dens.csv")
saveprettycsv(c2, pathmix, "mix_counts_eme.csv")
saveprettycsv(c3, pathmix, "mix_counts_emcp.csv")
saveprettycsv(c4, pathmix, "mix_counts_emv.csv")
saveprettycsv(c5, pathmix, "mix_counts_actcoeff.csv")
saveprettycsv(c6, pathmix, "mix_counts_sos.csv")
saveprettycsv(c7, pathmix, "mix_counts_dielec.csv")

# save pickle
saveprettypickle(dfbig, pathdf, "alldata_bin.pkl")
saveprettypickle(df1, pathdf, "dens_bin.pkl")
saveprettypickle(df2, pathdf, "eme_bin.pkl")
saveprettypickle(df3, pathdf, "emcp_bin.pkl")
saveprettypickle(df4, pathdf, "emv_bin.pkl")
saveprettypickle(df5, pathdf, "actcoeff_bin.pkl")
saveprettypickle(df6, pathdf, "sos_bin.pkl")
saveprettypickle(df7, pathdf, "dielec_bin.pkl")

saveprettypickle(a, pathjourn, "binname_counts_all.pkl")
saveprettypickle(a1, pathjourn, "binname_counts_dens.pkl")
saveprettypickle(a2, pathjourn, "binname_counts_eme.pkl")
saveprettypickle(a3, pathjourn, "binname_counts_emcp.pkl")
saveprettypickle(a4, pathjourn, "binname_counts_emv.pkl")
saveprettypickle(a5, pathjourn, "binname_counts_actcoeff.pkl")
saveprettypickle(a6, pathjourn, "binname_counts_sos.pkl")
saveprettypickle(a7, pathjourn, "binname_counts_dielec.pkl")

saveprettypickle(b, pathcomp, "bincomp_counts_all.pkl")
saveprettypickle(b1, pathcomp, "bincomp_counts_dens.pkl")
saveprettypickle(b2, pathcomp, "bincomp_counts_eme.pkl")
saveprettypickle(b3, pathcomp, "bincomp_counts_emcp.pkl")
saveprettypickle(b4, pathcomp, "bincomp_counts_emv.pkl")
saveprettypickle(b5, pathcomp, "bincomp_counts_actcoeff.pkl")
saveprettypickle(b6, pathcomp, "bincomp_counts_sos.pkl")
saveprettypickle(b7, pathcomp, "bincomp_counts_dielec.pkl")

saveprettypickle(c, pathmix, "mix_counts_all.pkl")
saveprettypickle(c1, pathmix, "mix_counts_dens.pkl")
saveprettypickle(c2, pathmix, "mix_counts_eme.pkl")
saveprettypickle(c3, pathmix, "mix_counts_emcp.pkl")
saveprettypickle(c4, pathmix, "mix_counts_emv.pkl")
saveprettypickle(c5, pathmix, "mix_counts_actcoeff.pkl")
saveprettypickle(c6, pathmix, "mix_counts_sos.pkl")
saveprettypickle(c7, pathmix, "mix_counts_dielec.pkl")


def SMILESchk(df, SMILES):
    ind = df.SMILES1.isin(SMILES) & df.SMILES2.isin(SMILES) 
    dfDavSet = df[ind]
    
    return dfDavSet

dfD = SMILESchk(df, S.SMILES)
df1D, a1D, b1D, c1D = dfpretty(dfD, "Mass density, kg/m3")
df2D, a2D, b2D, c2D = dfpretty(dfD, "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol")
df3D, a3D, b3D, c3D = dfpretty(dfD, "Excess molar heat capacity, J/K/mol")
df4D, a4D, b4D, c4D = dfpretty(dfD, "Excess molar volume, m3/mol")
df5D, a5D, b5D, c5D = dfpretty(dfD, "Activity coefficient")
df6D, a6D, b6D, c6D = dfpretty(dfD, "Speed of sound, m/s")
df7D, a7D, b7D, c7D = dfpretty(dfD, "Relative permittivity at zero frequency")

# save csv with ; delimiter
saveprettycsv(df1D, pathdf, "dens_bin_BDC.csv")
saveprettycsv(df2D, pathdf, "eme_bin_BDC.csv")
saveprettycsv(df3D, pathdf, "emcp_bin_BDC.csv")
saveprettycsv(df4D, pathdf, "emv_bin_BDC.csv")
saveprettycsv(df5D, pathdf, "actcoeff_bin_BDC.csv")
saveprettycsv(df6D, pathdf, "sos_bin_BDC.csv")
saveprettycsv(df7D, pathdf, "dielec_bin_BDC.csv")

saveprettycsv(a1D, pathjourn, "binname_counts_dens_BDC.csv")
saveprettycsv(a2D, pathjourn, "binname_counts_eme_BDC.csv")
saveprettycsv(a3D, pathjourn, "binname_counts_emcp_BDC.csv")
saveprettycsv(a4D, pathjourn, "binname_counts_emv_BDC.csv")
saveprettycsv(a5D, pathjourn, "binname_counts_actcoeff_BDC.csv")
saveprettycsv(a6D, pathjourn, "binname_counts_sos_BDC.csv")
saveprettycsv(a7D, pathjourn, "binname_counts_dielec_BDC.csv")

saveprettycsv(b1D, pathcomp, "comp_counts_dens_BDC.csv")
saveprettycsv(b2D, pathcomp, "comp_counts_eme_BDC.csv")
saveprettycsv(b3D, pathcomp, "comp_counts_emcp_BDC.csv")
saveprettycsv(b4D, pathcomp, "comp_counts_emv_BDC.csv")
saveprettycsv(b5D, pathcomp, "comp_counts_actcoeff_BDC.csv")
saveprettycsv(b6D, pathcomp, "comp_counts_sos_BDC.csv")
saveprettycsv(b7D, pathcomp, "comp_counts_dielec_BDC.csv")

saveprettycsv(c1D, pathmix, "mix_counts_dens_BDC.csv")
saveprettycsv(c2D, pathmix, "mix_counts_eme_BDC.csv")
saveprettycsv(c3D, pathmix, "mix_counts_emcp_BDC.csv")
saveprettycsv(c4D, pathmix, "mix_counts_emv_BDC.csv")
saveprettycsv(c5D, pathmix, "mix_counts_actcoeff_BDC.csv")
saveprettycsv(c6D, pathmix, "mix_counts_sos_BDC.csv")
saveprettycsv(c7D, pathmix, "mix_counts_dielec_BDC.csv")

# save pickle
saveprettypickle(df1D, pathdf, "dens_bin_BDC.pkl")
saveprettypickle(df2D, pathdf, "eme_bin_BDC.pkl")
saveprettypickle(df3D, pathdf, "emcp_bin_BDC.pkl")
saveprettypickle(df4D, pathdf, "emv_bin_BDC.pkl")
saveprettypickle(df5D, pathdf, "actcoeff_bin_BDC.pkl")
saveprettypickle(df6D, pathdf, "sos_bin_BDC.pkl")
saveprettypickle(df7D, pathdf, "dielec_bin_BDC.pkl")

saveprettypickle(a1D, pathjourn, "binname_counts_dens_BDC.pkl")
saveprettypickle(a2D, pathjourn, "binname_counts_eme_BDC.pkl")
saveprettypickle(a3D, pathjourn, "binname_counts_emcp_BDC.pkl")
saveprettypickle(a4D, pathjourn, "binname_counts_emv_BDC.pkl")
saveprettypickle(a5D, pathjourn, "binname_counts_actcoeff_BDC.pkl")
saveprettypickle(a6D, pathjourn, "binname_counts_sos_BDC.pkl")
saveprettypickle(a7D, pathjourn, "binname_counts_dielec_BDC.pkl")

saveprettypickle(b1D, pathcomp, "comp_counts_dens_BDC.pkl")
saveprettypickle(b2D, pathcomp, "comp_counts_eme_BDC.pkl")
saveprettypickle(b3D, pathcomp, "comp_counts_emcp_BDC.pkl")
saveprettypickle(b4D, pathcomp, "comp_counts_emv_BDC.pkl")
saveprettypickle(b5D, pathcomp, "comp_counts_actcoeff_BDC.pkl")
saveprettypickle(b6D, pathcomp, "comp_counts_sos_BDC.pkl")
saveprettypickle(b7D, pathcomp, "comp_counts_dielec_BDC.pkl")

saveprettypickle(c1D, pathmix, "mix_counts_dens_BDC.pkl")
saveprettypickle(c2D, pathmix, "mix_counts_eme_BDC.pkl")
saveprettypickle(c3D, pathmix, "mix_counts_emcp_BDC.pkl")
saveprettypickle(c4D, pathmix, "mix_counts_emv_BDC.pkl")
saveprettypickle(c5D, pathmix, "mix_counts_actcoeff_BDC.pkl")
saveprettypickle(c6D, pathmix, "mix_counts_sos_BDC.pkl")
saveprettypickle(c7D, pathmix, "mix_counts_dielec_BDC.pkl")






