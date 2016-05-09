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

df = th.pandas_dataframe()
dt = list(df.columns)

bad_filenames = ["/home/bmanubay/.thermoml/j.fluid.2013.12.014.xml"]  # This file confirmed to have possible data entry errors.
df = df[~df.filename.isin(bad_filenames)]

experiments = ["Mass density, kg/m3", "Partial molar enthalpy, kJ/mol", "Partial molar volume, m3/mol", "Partial pressure, kPa", "Excess molar Gibbs energy, kJ/mol", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol", "Excess molar heat capacity, J/K/mol", "Excess molar volume, m3/mol", "(Relative) activity", "Activity coefficient", "Speed of sound, m/s", "Relative permittivity at zero frequency"]

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

df["smiles1"] = df.x1.apply(lambda x: resolve_cached(x, "smiles"))  # This should be cached via sklearn.
df = df[df.smiles1 != None]
df = df[df["smiles1"].str.contains('=O') == False] # Getting rid of data sets with C=O and C=C occurrences
df = df[df["smiles1"].str.contains('#') == False]
df = df[df["smiles1"].str.contains('O=') == False]
df = df[df["smiles1"].str.contains('=C') == False]
df = df[df["smiles1"].str.contains('C=') == False]
df.dropna(subset=["smiles1"], inplace=True)
df = df.ix[df.smiles1.dropna().index]
df["smiles2"] = df.x2.apply(lambda x: resolve_cached(x, "smiles"))  # This should be cached via sklearn.
df = df[df.smiles2 != None]
df = df[df["smiles2"].str.contains('=O') == False] # Getting rid of data sets with C=O and C=C occurrences
df = df[df["smiles2"].str.contains('#') == False]
df = df[df["smiles2"].str.contains('O=') == False]
df = df[df["smiles2"].str.contains('=C') == False]
df = df[df["smiles2"].str.contains('C=') == False]
df.dropna(subset=["smiles2"], inplace=True)
df = df.ix[df.smiles2.dropna().index]

    
df["cas1"] = df.x1.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "cas")))  # This should be cached via sklearn.
df = df[df.cas1 != None]
df = df.ix[df.cas1.dropna().index]
df["cas2"] = df.x2.apply(lambda x: thermoml_lib.get_first_entry(resolve_cached(x, "cas")))  # This should be cached via sklearn.
df = df[df.cas2 != None]
df = df.ix[df.cas2.dropna().index]


# Neither names (components) nor smiles are unique.  Use CAS to ensure consistency.
cannonical_smiles_lookup1 = df.groupby("cas1").smiles1.first()
cannonical_components_lookup1 = df.groupby("cas1").x1.first()
cannonical_smiles_lookup2 = df.groupby("cas2").smiles2.first()
cannonical_components_lookup2 = df.groupby("cas2").x2.first()


df["smiles1"] = df.cas1.apply(lambda x: cannonical_smiles_lookup1[x])
df["x1"] = df.cas1.apply(lambda x: cannonical_components_lookup1[x])
df["smiles2"] = df.cas2.apply(lambda x: cannonical_smiles_lookup2[x])
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


dfbig = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"], df["Temperature, K"], df["Pressure, kPa"], df["Mass density, kg/m3"], df["Mass density, kg/m3_std"], df["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol"], df["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std"], df["Excess molar heat capacity, J/K/mol"], df["Excess molar heat capacity, J/K/mol_std"], df["Excess molar volume, m3/mol"], df["Excess molar volume, m3/mol_std"], df["Activity coefficient"], df["Activity coefficient_std"], df["Speed of sound, m/s"], df["Speed of sound, m/s_std"], df["Relative permittivity at zero frequency"], df["Relative permittivity at zero frequency_std"]], axis=1, keys=["filename", "x1", "x2", "smiles1", "smiles2", "components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Mass density, kg/m3", "Mass density, kg/m3_std", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std", "Excess molar heat capacity, J/K/mol", "Excess molar heat capacity, J/K/mol_std", "Excess molar volume, m3/mol", "Excess molar volume, m3/mol_std", "Activity coefficient", "Activity coefficient_std", "Speed of sound, m/s", "Speed of sound, m/s_std", "Relative permittivity at zero frequency", "Relative permittivity at zero frequency_std"])

dfbig.groupby(["filename"])

a = dfbig["filename"].value_counts()
a = a.reset_index()
a.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b = dfbig["x1"].value_counts()
b = b.reset_index()
b.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc = dfbig["smiles1"].value_counts()
cc = cc.reset_index()
cc.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b = pd.concat([b["Component"], cc["SMILES"], cc["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c = dfbig["x2"].value_counts()
c = c.reset_index()
c.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd = dfbig["smiles2"].value_counts()
dd = dd.reset_index()
dd.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c = pd.concat([c["Component"], dd["SMILES"], dd["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d = dfbig["components"].value_counts()
d = d.reset_index()
d.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df1 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Mass density, kg/m3"], df["Mass density, kg/m3_std"]], axis=1, keys=["filename", "x1", "x2", "smiles1", "smiles2","components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Mass density, kg/m3", "Mass density, kg/m3_std"])
df1["Mass density, kg/m3_std"].replace('nan', np.nan, inplace=True)
df1 = df1[np.isnan(df1["Mass density, kg/m3_std"])==False]
a1 = df1["filename"].value_counts()
a1 = a1.reset_index()
a1.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b1 = df1["x1"].value_counts()
b1 = b1.reset_index()
b1.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc1 = df1["smiles1"].value_counts()
cc1 = cc1.reset_index()
cc1.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b1 = pd.concat([b1["Component"], cc1["SMILES"], cc1["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c1 = df1["x2"].value_counts()
c1 = c1.reset_index()
c1.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd1 = df1["smiles2"].value_counts()
dd1 = dd1.reset_index()
dd1.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c1 = pd.concat([c1["Component"], dd1["SMILES"], dd1["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d1 = df1["components"].value_counts()
d1 = d1.reset_index()
d1.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df2 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol"], df["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std"]], axis=1, keys=["filename", "x1", "x2", "smiles1", "smiles2","components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol", "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std"])
df2["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std"].replace('nan', np.nan, inplace=True)
df2 = df2[np.isnan(df2["Excess molar enthalpy (molar enthalpy of mixing), kJ/mol_std"])==False]
a2 = df2["filename"].value_counts()
a2 = a2.reset_index()
a2.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b2 = df2["x1"].value_counts()
b2 = b2.reset_index()
b2.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc2 = df2["smiles1"].value_counts()
cc2 = cc2.reset_index()
cc2.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b2 = pd.concat([b2["Component"], cc2["SMILES"], cc2["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c2 = df2["x2"].value_counts()
c2 = c2.reset_index()
c2.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd2 = df2["smiles2"].value_counts()
dd2 = dd2.reset_index()
dd2.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c2 = pd.concat([c2["Component"], dd2["SMILES"], dd2["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d2 = df2["components"].value_counts()
d2 = d2.reset_index()
d2.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df3 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Excess molar heat capacity, J/K/mol"], df["Excess molar heat capacity, J/K/mol_std"]], axis=1, keys=["filename", "x1", "x2","smiles1", "smiles2", "components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Excess molar heat capacity, J/K/mol", "Excess molar heat capacity, J/K/mol_std"])
df3["Excess molar heat capacity, J/K/mol_std"].replace('nan', np.nan, inplace=True)
df3 = df3[np.isnan(df3["Excess molar heat capacity, J/K/mol_std"])==False]
a3 = df3["filename"].value_counts()
a3 = a3.reset_index()
a3.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b3 = df3["x1"].value_counts()
b3 = b3.reset_index()
b3.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc3 = df3["smiles1"].value_counts()
cc3 = cc3.reset_index()
cc3.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b3 = pd.concat([b3["Component"], cc3["SMILES"], cc3["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c3 = df3["x2"].value_counts()
c3 = c3.reset_index()
c3.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd3 = df3["smiles2"].value_counts()
dd3 = dd3.reset_index()
dd3.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c3 = pd.concat([c3["Component"], dd3["SMILES"], dd3["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d3 = df1["components"].value_counts()
d3 = d1.reset_index()
d3.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df4 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Excess molar volume, m3/mol"], df["Excess molar volume, m3/mol_std"]], axis=1, keys=["filename", "x1", "x2", "smiles1", "smiles2","components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Excess molar volume, m3/mol", "Excess molar volume, m3/mol_std"])
df4["Excess molar volume, m3/mol_std"].replace('nan', np.nan, inplace=True)
df4 = df4[np.isnan(df4["Excess molar volume, m3/mol_std"])==False]
a4 = df4["filename"].value_counts()
a4 = a4.reset_index()
a4.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b4 = df4["x1"].value_counts()
b4 = b4.reset_index()
b4.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc4 = df4["smiles1"].value_counts()
cc4 = cc4.reset_index()
cc4.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b4 = pd.concat([b4["Component"], cc4["SMILES"], cc4["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c4 = df1["x2"].value_counts()
c4 = c1.reset_index()
c4.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd4 = df4["smiles2"].value_counts()
dd4 = dd4.reset_index()
dd4.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c4 = pd.concat([c4["Component"], dd4["SMILES"], dd4["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d4 = df4["components"].value_counts()
d4 = d4.reset_index()
d4.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df5 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Activity coefficient"], df["Activity coefficient_std"]], axis=1, keys=["filename", "x1", "x2", "smiles1", "smiles2","components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Activity coefficient", "Activity coefficient_std"])
df5["Activity coefficient_std"].replace('nan', np.nan, inplace=True)
df5 = df5[np.isnan(df5["Activity coefficient_std"])==False]
a5 = df5["filename"].value_counts()
a5 = a5.reset_index()
a5.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b5 = df5["x1"].value_counts()
b5 = b5.reset_index()
b5.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc5 = df5["smiles1"].value_counts()
cc5 = cc5.reset_index()
cc5.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b5 = pd.concat([b5["Component"], cc5["SMILES"], cc5["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c5 = df1["x2"].value_counts()
c5 = c5.reset_index()
c5.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd5 = df5["smiles2"].value_counts()
dd5 = dd5.reset_index()
dd5.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c5 = pd.concat([c5["Component"], dd5["SMILES"], dd5["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d5 = df5["components"].value_counts()
d5 = d5.reset_index()
d5.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df6 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"], df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Speed of sound, m/s"], df["Speed of sound, m/s_std"]], axis=1, keys=["filename", "x1", "x2","smiles1", "smiles2", "components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Speed of sound, m/s", "Speed of sound, m/s_std"])
df6["Speed of sound, m/s_std"].replace('nan', np.nan, inplace=True)
df6 = df6[np.isnan(df6["Speed of sound, m/s_std"])==False]
a6 = df6["filename"].value_counts()
a6 = a6.reset_index()
a6.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b6 = df6["x1"].value_counts()
b6 = b6.reset_index()
b6.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc6 = df6["smiles1"].value_counts()
cc6 = cc6.reset_index()
cc6.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b6 = pd.concat([b6["Component"], cc6["SMILES"], cc6["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c6 = df6["x2"].value_counts()
c6 = c6.reset_index()
c6.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd6 = df6["smiles2"].value_counts()
dd6 = dd6.reset_index()
dd6.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c6 = pd.concat([c6["Component"], dd6["SMILES"], dd6["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d6 = df6["components"].value_counts()
d6 = d6.reset_index()
d6.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

df7 = pd.concat([df['filename'], df['x1'], df['x2'], df["smiles1"], df["smiles2"],df["components"], df["Mole fraction"],df["Temperature, K"], df["Pressure, kPa"], df["Relative permittivity at zero frequency"], df["Relative permittivity at zero frequency_std"]], axis=1, keys=["filename", "x1", "x2","smiles1", "smiles2", "components", "Mole fraction", "Temperature, K", "Pressure, kPa", "Relative permittivity at zero frequency", "Relative permittivity at zero frequency_std"])
df7["Relative permittivity at zero frequency_std"].replace('nan', np.nan, inplace=True)
df7 = df7[np.isnan(df7["Relative permittivity at zero frequency_std"])==False]
a7 = df7["filename"].value_counts()
a6 = a6.reset_index()
a6.rename(columns={"index":"Filename","filename":"Count"},inplace=True)
b7 = df7["x1"].value_counts()
b7 = b7.reset_index()
b7.rename(columns={"index":"Component","x1":"Count"},inplace=True)
cc7 = df7["smiles1"].value_counts()
cc7 = cc7.reset_index()
cc7.rename(columns={"index":"SMILES","smiles1":"Count"},inplace=True)
b7 = pd.concat([b7["Component"], cc7["SMILES"], cc7["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
c7 = df7["x2"].value_counts()
c7 = c7.reset_index()
c7.rename(columns={"index":"Component","x2":"Count"},inplace=True)
dd7 = df7["smiles2"].value_counts()
dd7 = dd7.reset_index()
dd7.rename(columns={"index":"SMILES","smiles2":"Count"},inplace=True)
c7 = pd.concat([c7["Component"], dd7["SMILES"], dd7["Count"]], axis=1, keys=["Component", "SMILES", "Count"])
d7 = df7["components"].value_counts()
d7 = d7.reset_index()
d7.rename(columns={"index":"Mixture","components":"Count"},inplace=True)

dfbig.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/Ken_binary_sets_all.csv")
a.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_all.csv")
b.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_all.csv")
c.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_all.csv")
d.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_all.csv")

df1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/dens_bin.csv")
a1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_dens.csv")
b1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_dens.csv")
c1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_dens.csv")
d1.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_dens.csv")

df2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/eme_bin.csv")
a2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_eme.csv")
b2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_eme.csv")
c2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_eme.csv")
d2.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_eme.csv")

df3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/emcp_bin.csv")
a3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_emcp.csv")
b3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_emcp.csv")
c3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_emcp.csv")
d3.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_emcp.csv")

df4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/emv_bin.csv")
a4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_emv.csv")
b4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_emv.csv")
c4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_emv.csv")
d4.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_emv.csv")

df5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/actcoeff_bin.csv")
a5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_act.csv")
b5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_act.csv")
c5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_act.csv")
d5.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_act.csv")

df6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/sos_bin.csv")
a6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_sos.csv")
b6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_sos.csv")
c6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_sos.csv")
d6.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_sos.csv")

df7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Property data/dielec_bin.csv")
a7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Journal name counts/name_counts_dielec.csv")
b7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x1_counts_dielec.csv")
c7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Component counts/x2_counts_dielec.csv")
d7.to_csv("/home/bmanubay/.thermoml/tables/Ken/Binary/Mixture counts/mix_counts_dielec.csv")
