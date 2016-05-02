
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
import openeye.oedepict as oedepict
from IPython.display import display


# In[3]:

def depict(mol, width=500, height=200):
    from IPython.display import Image
    dopt = oedepict.OEPrepareDepictionOptions()
    dopt.SetDepictOrientation( oedepict.OEDepictOrientation_Horizontal)
    oedepict.OEPrepareDepiction(mol, dopt)
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    ofs = oechem.oeosstream()
    oedepict.OERenderMolecule(ofs, 'png', disp)
    ofs.flush()
    return Image(data = "".join(ofs.str()))


# In[15]:

# filters for peroxides, 3-memb-rings, ketals
badSmarts = [ 'OO', '[r3]', '[O&H1]C[O&H1]' ]
badSSlst = [ oechem.OESubSearch( smarts) for smarts in badSmarts]


# In[16]:

fileprefix= 'AlkEthOH_chain'
#fileprefix= 'AlkEthOH_rings'
ifs = oechem.oemolistream(fileprefix+'.smi')
ofs = oechem.oemolostream(fileprefix+'_filt1.smi')
#molLst1 = [ mol for mol in ifs.GetOEMols() ]


# In[17]:

mol = oechem.OEMol()
for mol in ifs.GetOEMols():
    goodMol = True
    for badSS in badSSlst:
        oechem.OEPrepareSearch(mol, badSS)
        if badSS.SingleMatch(mol):
            goodMol = False
            display( depict(mol))
    if goodMol:
        oechem.OEWriteConstMolecule(ofs, mol)


# In[18]:

ifs.close()
ofs.close()


# In[ ]:



