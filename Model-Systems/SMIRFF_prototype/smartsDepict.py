
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

def depictMatch(mol, match, width=500, height=200):
    from IPython.display import Image
    dopt = oedepict.OEPrepareDepictionOptions()
    dopt.SetDepictOrientation( oedepict.OEDepictOrientation_Horizontal)
    dopt.SetSuppressHydrogens(True)
    oedepict.OEPrepareDepiction(mol, dopt)
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    hstyle = oedepict.OEHighlightStyle_Color
    hcolor = oechem.OEColor(oechem.OELightBlue)
    oedepict.OEAddHighlighting(disp, hcolor, hstyle, match)
    ofs = oechem.oeosstream()
    oedepict.OERenderMolecule(ofs, 'png', disp)
    ofs.flush()
    return Image(data = "".join(ofs.str()))


# In[4]:

Smarts = '[#6X4]-[#6X4]-[#8X2]'
qmol = oechem.OEQMol()
if not oechem.OEParseSmarts( qmol, Smarts ):
    print( 'OEParseSmarts failed')
ss = oechem.OESubSearch( qmol)


# In[5]:

fileprefix= 'AlkEthOH_dvrs1'
ifs = oechem.oemolistream(fileprefix+'.oeb')


# In[6]:

mol = oechem.OEMol()
for mol in ifs.GetOEMols():
    goodMol = True
    oechem.OEPrepareSearch(mol, ss)
    unique = True
    for match in ss.Match(mol, unique):
        display( depictMatch(mol, match))
        for ma in match.GetAtoms():
            print(ma.target.GetIdx(), end=" ")
            #print(ma.pattern.GetIdx(), end=" ")


# In[7]:

ifs.close()


# In[7]:



