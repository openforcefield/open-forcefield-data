# open-forcefield-data/Initial Molecule Choices
Organized csv files showing data point count statistics of molecules and mixtures across all pure solvent and binary mixture properties

All files with '_all' or lone '_interesting' title appendage were created using the thermomlcnts.py script in the 'Python' directory (https://github.com/open-forcefield-group/open-forcefield-data/tree/master/Python)

Any files with the '_diverse' title appendage (or any combination of it) was created using the diversitychk.py script in the 'Python' directory (https://github.com/open-forcefield-group/open-forcefield-data/tree/master/Python)

# MAIN csv contents ('_all' title appendage) 

allcomp_counts_all.csv - counts of number of data points for all molecules appearing in the pure and binary mixture sets across all properties (no additional filters applied)

bincomp_counts_all.csv - Count of data points per unique molecule across all binary properties. It’s important to note that while this count was performed on individual components, it is for the binary property sets (no additional filters applied)

mix_counts_all.csv - count of the number of data points per unique binary mixture across all mixture properties (no additional filters applied)

purecomp_counts_all.csv - Count of data points per unique molecule across all pure properties (no additional filters applied)

# Additional title appendange meanings

_diverse - These data point counts represent only molecules on Chris's diverse list (see: https://github.com/open-forcefield-group/open-forcefield-data/blob/master/Model-Systems/AlkEthOH_distrib/AlkEthOH_chain_filt1.smi and https://github.com/open-forcefield-group/open-forcefield-data/blob/master/Model-Systems/AlkEthOH_distrib/AlkEthOH_rings_filt1.smi) 

_interesting - These data point counts contain only binary data for those mixtures which are not alkane-alkane

**Both appendages would mean both filters were applied
