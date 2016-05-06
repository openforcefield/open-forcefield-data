
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

#fileprefix= 'AlkEthOH_rings_filt1'
#fileprefix= 'test_filt1'
fileprefix= 'AlkEthOH_chain_filt1'
ifs = oechem.oemolistream(fileprefix+'.smi')
ofs = oechem.oemolostream(fileprefix+'_namedCharged.oeb')


# In[4]:

maxConfs = 1
omega = oeomega.OEOmega()
omega.SetMaxConfs( maxConfs)
omega.SetStrictStereo( False)
mol = oechem.OEMol()
for mol in ifs.GetOEMols():
    #print( mol.GetTitle() )
    if not omega(mol):
        print( "omega failed on %s" % mol.GetTitle() )
    title= mol.GetTitle()
    oechem.OETriposAtomNames( mol)
    oequacpac.OEAssignPartialCharges(mol, oequacpac.OECharges_AM1BCCSym)
    oechem.OEWriteConstMolecule(ofs, mol)


# In[5]:

ifs.close()
ofs.close()


# In[6]:

xpatParmaFrosst = '/Users/bayly/BaylyData/collaborations/sabbatical/parmaFrosst/parmaFrosst_atyper_v1.2.xpat'
cmd_string = 'xpat -p '+xpatParmaFrosst+' -in '+fileprefix+'_namedCharged.oeb'+' -out '+fileprefix+'.oeb'
print( cmd_string)
os.system( cmd_string)


# In[ ]:



