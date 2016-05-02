
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
get_ipython().magic(u'matplotlib inline')


# In[2]:

import openeye.oechem as oechem
import openeye.oeszybki as oeszybki
import openeye.oequacpac as oequacpac
import openeye.oeomega as oeomega
import openeye.oedepict as oedepict
from IPython.display import display


# In[3]:

def hasAmberParams( mol, cmd_string):
    ofslig = oechem.oemolostream( 'lig.mol2')
    ofslig.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield );
    oechem.OEWriteConstMolecule(ofslig, mol)
    ofslig.close()
    os.system( cmd_string)
    paramsNotSaved = 'Parameter file was not saved'
    leaplog = open( 'leap_lig.stdout', 'r' ).read()
    hasParams = not paramsNotSaved in leaplog
    return hasParams


# In[4]:

fileprefix= 'AlkEthOH_chain_filt1'
#fileprefix= 'AlkEthOH_rings_filt1'
#fileprefix= 'test_filt1'
#fileprefix= 'AlkEthOH_r21'
#fileprefix= 'thp23diol'
ifs = oechem.oemolistream(fileprefix+'.oeb')
ofs = oechem.oemolostream(fileprefix+'_hasParam.oeb')
ofsFail = oechem.oemolostream(fileprefix+'_lacksParam.oeb')
#ofs.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield );


# In[5]:

cmd_string = 'tleap -f leaprc.Frosst_AlkEthOH -f lig.leap_in >| leap_lig.stdout'
mol = oechem.OEMol()
for mol in ifs.GetOEMols():
    if hasAmberParams(mol,cmd_string):
        print( 'Has   Amber params: %s' % mol.GetTitle() )
        oechem.OEWriteConstMolecule(ofs, mol)
    else:
        print( 'Lacks Amber params: %s' % mol.GetTitle() )
        oechem.OEWriteConstMolecule(ofsFail, mol)


# In[6]:

ifs.close()
ofs.close()
ofsFail.close()


# In[ ]:



