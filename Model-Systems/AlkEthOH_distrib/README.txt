Sun May  1 17:51:15 MDT 2016
author: Christopher I. Bayly

This README.txt file describes the contents of this directory and
how to use them.

There are 3 key control files for use with amber:
-rw-r--r--  1 bayly  admin    4866 Apr 30 23:25 frcmod.Frosst_AlkEthOH
-rw-r--r--  1 bayly  admin     664 Apr 30 23:15 leaprc.Frosst_AlkEthOH
-rw-r--r--  1 bayly  staff    5292 Mar 12  2013 parmaFrosst_atyper_v1.2.xpat
The leaprc file leaprc.Frosst_AlkEthOH needs to be copied into
the amber tree directory $AMBERHOME/dat/leap/cmd/ and the
frcmod file frcmod.Frosst_AlkEthOH needs to be copied into
the amber tree directory $AMBERHOME/dat/leap/parm/. The xpat file
parmaFrosst_atyper_v1.2.xpat is given here for reference purposes
to show the SMARTS strings used for parm@Frosst atom typing,
unfortunately the xpat executable cannot be distributed at this time.

There are 4 python scripts forming an overall workflow which will be
described below:
-rw-r--r--@ 1 bayly  staff    1959 May  1 15:37 oebHasAmberParams.py
-rw-r--r--@ 1 bayly  staff    2097 May  1 17:47 oebToAmberMol2TopCrd.py
-rw-r--r--@ 1 bayly  staff    1629 May  1 15:39 smiAmberize.py
-rw-r--r--@ 1 bayly  staff    1815 May  1 15:03 substrucRemove.py
These were dumped straight out of iPython notebooks so they may not
run standalone as is, but they clearly show the workflow I used.

There are 3 molecule sets beginning as smiles, with pdfs of the structures:
-rw-r--r--  1 bayly  staff   38894 Apr 30 16:35 AlkEthOH_chain.smi
-rw-r--r--  1 bayly  staff  285637 Apr 30 16:38 AlkEthOH_chain.pdf
-rw-r--r--  1 bayly  staff   34105 Apr 30 16:34 AlkEthOH_rings.smi
-rw-r--r--  1 bayly  staff  292609 Apr 30 16:38 AlkEthOH_rings.pdf
-rw-r--r--  1 bayly  staff    1311 Apr 30 16:58 test.smi
-rw-r--r--@ 1 bayly  staff   12401 May  1 19:37 test.pdf
The files AlkEthOH_chain.* are acyclic alkanes, alcohols, and ethers
up to about 8 or so carbons. The files AlkEthOH_rings.* are cyclic
alkanes, alcohols, and ethers with rings sizes varying from 3 to 9.
The python script used to generate these is not included, but while
intended to generate a large and diverse set it is far from exhaustive.
The goal was to sample small rings, sterically congested substitutions,
and dense functionalisation so as to exercise the parameterization
scheme, the non-bonded parameters, and the torsion parameters. File
test.smi has the top 50 smiles from AlkEthOH_rings.smi, generated
for testing purposes... a good small set to prototype a workflow.

There are 3 final oeb files ready to be converted into amber .top and .crd
files, with pdfs of the structures:
-rw-r--r--  1 bayly  staff  661638 May  1 15:26 AlkEthOH_chain_filt1.oeb
-rw-r--r--  1 bayly  staff  758310 Apr 30 20:18 AlkEthOH_rings_filt1.oeb
-rw-r--r--  1 bayly  staff   24464 Apr 30 20:08 test_filt1.oeb
These are the work product of the workflow detailed below, ready to
be converted into amber input files, having one fully-protonated geometry
for each structure, with atom names, AM1BCC charges, and parm@Frosst
atom types from xpat.

Finally, there are several intermediate files generated along the way:
-rw-r--r--@ 1 bayly  staff  186884 May  1 15:43 AlkEthOH_chain_filt1.pdf
-rw-r--r--  1 bayly  staff   26904 May  1 15:02 AlkEthOH_chain_filt1.smi
-rw-r--r--@ 1 bayly  staff  264271 Apr 30 17:52 AlkEthOH_rings_filt1.pdf
-rw-r--r--  1 bayly  staff   30967 Apr 30 17:51 AlkEthOH_rings_filt1.smi
-rw-r--r--@ 1 bayly  staff   11287 May  1 15:52 test_filt1.pdf
-rw-r--r--  1 bayly  staff    1128 Apr 30 17:50 test_filt1.smi


Additional files added by Victoria T. Lim on 24 June 2016
The directory files_for_c1302/ contains 2 files particular to 
generating AMBER input files for AlkEthOH_c1302 (water).
(1) AlkEthOH_c1302_edited.leap_in
    tleap input file with extra line to read in frcmod.tip3p
(2) frcmod.tip3p
    parameters for tip3p water, not included in frcmod.Frosst_AlkEthOH


Workflow
--------
To demonstrate the workflow I will start with the small test set
-rw-r--r--  1 bayly  staff    1311 Apr 30 16:58 test.smi

1. Remove undesired substructures.
Python script substrucRemove.py is run first to identify structures
containing peroxides (-O-O-), 3-membered rings, and ketals (HO-Csp3-OH),
all of which cause too much complication for now. Substituting 'test'
for 'AlkEthOH_chain' in the following line of the script:
fileprefix= 'AlkEthOH_chain'
will read in test.smi and generate file
-rw-r--r--  1 bayly  staff    1128 Apr 30 17:50 test_filt1.smi
which has the offending chemistry removed. This file is carried forward.

2. Amberize with atom names, charges, and parm@Frosst atom types.
Python script smiAmberize.py take a smiles string, generates a single
fully-protonated conformer, names the atoms with Tripos atom names
(fine for our purposes, any will do to satisfy .mol2 format),
charges it with single-conformer AM1BCC charges (adequate for now),
and generates parm@Frosst atom types for it by running xpat as a system command.
Substituting 'test_filt1' for 'AlkEthOH_chain_filt1' in the following
line of the script:
fileprefix= 'AlkEthOH_chain_filt1'
will read in test.smi and generate file:
-rw-r--r--  1 bayly  staff   24464 Apr 30 20:08 test_filt1.oeb
It will also generate an intermediate file test_filt1_namedCharged.oeb
as the input to the xpat system command, but that file was not included
here.

3. Check that there are no missing amber parameters in the frcmod file.
Python script oebHasAmberParams.py goes through a file of many
structures just to check if there are missing amber parameters in the
frcmod file you wish to parameterize it with. 
Substituting 'test_filt1' for 'AlkEthOH_chain_filt1' in the following
line of the script:
fileprefix= 'AlkEthOH_chain_filt1'
will generate two files, test_filt1_hasParam.oeb and
test_filt1_lacksParam.oeb, neither of which were included here because
in this case test_filt1_lacksParam.oeb was empty and test_filt1_hasParam.oeb
was identical to the input file because all structures have parameters.

4. Generate the amber input .mol2, .top, and .crd files
Python script oebToAmberMol2TopCrd.py generates, for each molecule
individually, the input .mol2 needed for tleap and the tleap input file
itself. It then runs tleap as a system command, using the
leaprc.Frosst_AlkEthOH and frcmod.Frosst_AlkEthOH specified in the
system command, to generate the amber-ready .top and .crd files.
The 'test_filt1' filename prefix is already correctly specified in
this script; it generates 168 files: 42 molecules times 4 files per
molecule (.mol2, .leap_in, .top, and .crd).
NOTE WELL that this could generate a large number of files (~4000)
in the directory for AlkEthOH_chain_filt1.oeb or AlkEthOH_rings_filt1.oeb.

* The oebToAmberMol2TopCrd.py was modified by VTL to run on all AlkEthOH
sets (chain, rings, test). The files are placed in a new dir (inputfiles/)
with subdirectories for each set. The number of files generated are:
chain (913*4 = 3652), ring (1036*4 = 4144), test (42*4 = 168).

This workflow was run on AlkEthOH_chain.smi and AlkEthOH_rings.smi
to produce amber-ready files:
-rw-r--r--  1 bayly  staff  661638 May  1 15:26 AlkEthOH_chain_filt1.oeb
-rw-r--r--  1 bayly  staff  758310 Apr 30 20:18 AlkEthOH_rings_filt1.oeb
which were checked in step 3 of the workflow to function with
the parm@Frosst parameter list abridged for alkanes, alcohols, and ethers.


parm@Frosst for alkanes, alcohols, and ethers
---------------------------------------------
The file frcmod.Frosst_AlkEthOH is the alkane, alcohol, and ether
subspace of the parm99/parm@Frosst parameter set, with the exception
of 3-membered rings (needing carbon atom type CJ) and ketals
(needing hydrogen atom type CX). This subspace is suitable for
toy dataset purposes, comprising very few atom types, being (from
the leaprc.Frosst_AlkEthOH file):
logFile leap_lig.log
#
addAtomTypes {
  { "HO"  "H" "sp3" }
  { "H1"  "H" "sp3" }
  { "H2"  "H" "sp3" }
  { "H3"  "H" "sp3" }
  { "HC"  "H" "sp3" }
  { "OH"  "O" "sp3" }
  { "OS"  "O" "sp3" }
  { "CT"  "C" "sp3" }
}

These few atom types make for a tractable number of parameters to
work with in the toy dataset. Tractable here is formally 72
groups, where a group corresponds to all the parameters needed for
a certain degree of freedom: 2 for bonds, 2 for angles, 3 per
trig function for torsions, 2 for vdW, and 1 for mass.


This concludes this README file. Please contact me with questions.
