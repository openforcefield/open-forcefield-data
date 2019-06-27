# Examples from the PDB database (https://www.rcsb.org/)
Ligands were extracted using the inhouse Proasis
software and filtered by multiple critiriea:

- Electron density quality criteria such es Resolution and Rvalue and others
- MW <= 500, Rotatable Bonds <= 6, no none-organic atoms

## pubLigsNeutralGoodDensity.sdf.gz:
648 ligands from pdb. Starting with criteriea mentioend above then further processed as follows:

- Computed pKa: most basic > 6, most acidic < 8
- All ionized atoms where neutralized by adding or removing H when possible. 
- Some cleanup of the connectivity was performed but there might still be issues.


## pubLigsChargedGoodDensity.sdf.gz
382 ligands from pdb. Starting with criteriea mentioend above then further processed as follows:

- Computed pKa: most basic > 8, most acidic < 6
- Ligands where also computationally protonated to their assumed state at pH 7.0. 
  Note: the computational method used is known to generate wrong protonation states and tautomers with
low population.
- Some cleanup of the connectivity was performed but there might still be issues.
