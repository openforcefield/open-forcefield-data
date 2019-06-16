# Utilize All Parameters

This is an attempt to create a set of (preferably) small molecules which utilize
all parameters in the smirnoff99Frosst force field. Currently, the
`smirnoff_utilize_all_params.ipynb` notebook shows several attempts which
utilize a set of molecules from Roche, as well as the first 1,000 and 10,000
molecules (these datasets are all included). However, none of these attempts
succeeded in covering the entire set of parameters, each leaving out over 100
parameters.

Hence, we decided to convert the notebook to a script,
`check_param_coverage.py`, which can be run on a computing cluster such as Green
Planet to produce parameter ids for a much larger set of molecules. We have done
one run so far over a million molecules split across 10 jobs, but 6 of those
failed and will be re-run.

Next, we created `select_molecules.py` to apply a greedy weighted set cover (see
[here](https://www.cs.huji.ac.il/course/2005/algo2/scribes/lecture2.pdf) for
more details) to find molecules that are small but cover as many parameters as
possible.

Finally, we calculated parameters for 1 million molecules in eMolecules and
selected a set of molecules from it. We ended up with a set of 68 molecules
which cover 305 parameters in smirnoff. The results are stored in the `selected`
folder, which contains the following files:
- `chosen.smi` - the SMILES strings of the chosen molecules
- `param_ids.json` - the list of parameter IDs of the chosen molecules
- `remaining_ids.json` - the list of remaining parameter IDs
- `params_by_molecule.json` - a mapping from the index of each molecule in
  `chosen.smi` to the molecule's SMILES string and the parameter IDs in that
  molecule.  Each ID is annotated with the list of atom indices where the
  parameter appeared in the molecule (there are multiple lists if the parameter
  appeared multiple times).

## Running

Before running any scripts here, set up the conda environment with
```
conda env create -f environment.yml
```
The scripts have the following dependencies:
- openeye-toolkits 2019.4.2
- openforcefield 0.4.0
