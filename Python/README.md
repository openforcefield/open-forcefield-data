# open-forcefield-data/Python
Python scripts for organizing and parsing ThermoML data using ThermoPyL

# Script contents and use

Ken_test_score.py - Script used to rank journal articles to submit to Ken for red flagging and uncertainty estimation. Scoring rewards large number of components per journal and wide distribution of data across properties measured.

binsrch.py - Script used to generate binary mixture property sets for all potential molecules/mixtures

diversitychk.py - Script used to check which of our initial molecules (see 'allcomp_counts_all.csv' in 'Initial Molecule Choices' directory) are in Chris's diverse set and generate allcomp_counts_diverse.csv (also in 'Initial Molecule Choices' directory)

puresrch.py - Script used to generate pure solvent property sets for all potential molecules

thermomlcnts.py - Script used to generate data point count data for all property sets passing initial atom type filter (See contents of 'Initial Molecule Choices' for outputs)



