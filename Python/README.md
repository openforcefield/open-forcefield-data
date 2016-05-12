# open-forcefield-data/Python
Python scripts for organizing and parsing ThermoML data using ThermoPyL

# Script contents and use

Asetcounts.py - Script used to generate data point counts for initial A set (Now degenerate due to change of initial ideas for A set. Still useful for showing ThermoPyL functionality)

binsrch.py - Script used to generate binary mixture property sets for all potential molecules/mixtures

diversitychk.py - Script used to check which of our initial molecules (see 'allcomp_counts_all.csv' in 'Initial Molecule Choices' directory) are in Chris's diverse set and generate allcomp_counts_diverse.csv (also in 'Initial Molecule Choices' directory)

mixchk.py - Script used to eliminate alkane-alkane mixtures from the 'mix_counts_all.csv' statistics. List is also prefiltered to eliminate any mixture containing a molecule not in Chris's diverse list. Output is 'mix_counts_interesting.csv' in 'Initial Molecule Choices' directory.

puresrch.py - Script used to generate pure solvent property sets for all potential molecules

thermomlcnts.py - Script used to generate data point count data for all property sets passing initial atom type filter (See contents of 'Initial Molecule Choices' for outputs)



