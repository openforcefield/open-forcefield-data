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

## Next Steps

We are still missing parameters, so we plan to use substructure searches on the
eMolecules database to find molecules which match the parameters we are missing.
We will also finish gathering parameters from more molecules in eMolecules.

## Running

Before running any scripts here, set up the conda environment with
```
conda env create -f environment.yml
```
The scripts have the following dependencies:
- openeye-toolkits 2019.4.2
- openforcefield 0.4.0


