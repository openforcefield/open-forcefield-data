
# coding: utf-8

# In[1]:

#!/usr/bin/env python
################################################################
#  Copyright (C) 2015 OpenEye Scientific Software, Inc.
################################################################

from __future__ import print_function
import os, sys
import pandas as pd
import scipy.stats as stats
import scipy as sci
import numpy as np
import pylab


# In[2]:

import openeye.oechem as oechem
import openeye.oeszybki as oeszybki
import openeye.oequacpac as oequacpac
import openeye.oeomega as oeomega
import openeye.oedepict as oedepict
from IPython.display import display


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

fileprefix= 'AlkEthOH_rings_filt1'
#fileprefix= 'AlkEthOH_r21'
ifs = oechem.oemolistream(fileprefix+'.oeb')


# In[5]:

cmd_string = 'leaprc.Frosst_AlkEthOH'
mol = oechem.OEMol()
for mol in ifs.GetOEMols():
    if hasAmberParams(mol,cmd_string):
        print( '%s successful writing amber .mol2, .top, and .crd file' % mol.GetTitle() )
    else:
        print( '%s failed writing amber .top file' % mol.GetTitle() )


# In[6]:

ifs.close()


# In[ ]:



