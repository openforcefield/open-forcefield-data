
# coding: utf-8

# In[1]:

#!/usr/bin/env python
################################################################
#  Copyright (C) 2015 OpenEye Scientific Software, Inc.
################################################################

# run this script in AlkEthOH_distrib directory containing *.oeb

from __future__ import print_function
import os, sys
from shutil import copyfile


# In[2]:

import openeye.oechem as oechem
import openeye.oeszybki as oeszybki
import openeye.oequacpac as oequacpac
import openeye.oeomega as oeomega
import openeye.oedepict as oedepict
#from IPython.display import display


# In[3]:

def hasAmberParams( mol, leaprc_string):
    # write mol2 file
    title = mol.GetTitle()
    ofsmol2 = oechem.oemolostream(title+'.mol2')
    ofsmol2.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield );
    oechem.OEWriteConstMolecule(ofsmol2, mol)
    ofsmol2.close()
    # write tleap input file
    os.system( 'echo """lig = loadMol2 '+title+'.mol2""" >| '+title+'.leap_in')
    os.system( 'echo """saveAmberParm lig '+title+'.top '+title+'.crd""" >> '+title+'.leap_in')
    os.system( 'echo quit >> '+title+'.leap_in')
    # run tleap
    #print( 'tleap -f '+leaprc_string+' -f '+title+'.leap_in >| leap_lig.stdout' )
    os.system( 'tleap -f '+leaprc_string+' -f '+title+'.leap_in >| leap_lig.stdout' )
    # check if param file was not saved (implies parameterization problems)
    paramsNotSaved = 'Parameter file was not saved'
    leaplog = open( 'leap_lig.stdout', 'r' ).read()
    hasParams = not paramsNotSaved in leaplog
    return hasParams


# In[4]:

fileprefixes = ['AlkEthOH_chain_filt1','AlkEthOH_rings_filt1','test_filt1']
cmd_string = 'leaprc.Frosst_AlkEthOH'

# create directory for all input files
if not os.path.exists('inputfiles'):
    os.makedirs('inputfiles')
os.chdir('inputfiles')



for fileprefix in fileprefixes: # run for all three .oeb files

    # create subdirectory for this set
    if not os.path.exists(fileprefix):
        os.makedirs(fileprefix)
    os.chdir(fileprefix)

    # copy temporary files
    copyfile('../../frcmod.Frosst_AlkEthOH','./frcmod.Frosst_AlkEthOH')
    copyfile('../../leaprc.Frosst_AlkEthOH','./leaprc.Frosst_AlkEthOH')
    copyfile('../../'+fileprefix+'.oeb','./'+fileprefix+'.oeb')


    ifs = oechem.oemolistream(fileprefix+'.oeb')
    mol = oechem.OEMol()
    for mol in ifs.GetOEMols():
        # add atom names c0 (methane) and c1302 (water) 
        if (oechem.OECount(mol, oechem.OEIsHeavy()) == 1):
            oechem.OETriposAtomNames(mol)
        # generate input files
        if hasAmberParams(mol,cmd_string):
            print( '%s successful writing amber .mol2, .top, and .crd file' % mol.GetTitle() )
            
        # treat water with diff pre-existing tleap input file
        elif mol.GetTitle().split("_")[1] == 'c1302':
            copyfile('../../files_for_c1302/frcmod.tip3p','./frcmod.tip3p')
            copyfile('../../files_for_c1302/AlkEthOH_c1302_edited.leap_in','./AlkEthOH_c1302_edited.leap_in')           
            os.system( 'tleap -f leaprc.Frosst_AlkEthOH -f AlkEthOH_c1302_edited.leap_in >| leap_lig.stdout' )
            print( '%s successful writing amber .mol2, .top, and .crd file' % mol.GetTitle() )

            os.remove('frcmod.tip3p')
            os.remove('AlkEthOH_c1302_edited.leap_in')

        else:                                                                      
            print( '%s failed writing amber .top file' % mol.GetTitle() ) 

    ifs.close()
    os.remove('frcmod.Frosst_AlkEthOH')
    os.remove('leaprc.Frosst_AlkEthOH')
    os.remove(fileprefix+'.oeb')
    os.chdir('../')





