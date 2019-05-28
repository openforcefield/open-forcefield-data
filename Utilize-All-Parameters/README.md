# Utilize All Parameters

This is an attempt to create a set of (preferably) small molecules which utilize
all parameters in the smirnoff99Frosst force field. Currently, the
`smirnoff_utilize_all_params.ipynb` notebook shows several attempts which
utilize a set of molecules from Roche, as well as the first 1,000 and 10,000
molecules (these datasets are all included). However, none of these attempts
succeeded in covering the entire set of parameters, each leaving out over 100
parameters.

## Next Steps

We will convert the notebook to a script and run it on 1 random sample of 1
million molecules from eMolecules. This will allow us to have a large set that
(hopefully) contains all parameters from smirnoff99Frosst.
