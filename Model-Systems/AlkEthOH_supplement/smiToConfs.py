#!/usr/bin/env python
#  generate conformers from smiles and assign AM1BCC charges and atom names
################################################################
import os, sys
import openeye.oechem as oechem
import openeye.oequacpac as oequacpac
import openeye.oeomega as oeomega

def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)

    if oechem.OECheckHelp(itf, argv):
        return 0

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    ifs = oechem.oemolistream()
    if not ifs.open(itf.GetString("-i")):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % itf.GetString("-i"))

    ofs = oechem.oemolostream()
    if not ofs.open(itf.GetString("-o")):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % itf.GetString("-o"))
    
    omega = oeomega.OEOmega()
    maxConfs = 0
    omega.SetMaxConfs( maxConfs)
    omega.SetStrictStereo( False)
    omega.SetSampleHydrogens( True)
    omega.SetEnumNitrogen( oeomega.OENitrogenEnumeration_All)
    mol = oechem.OEMol()
    for mol in ifs.GetOEMols():
        #print( mol.GetTitle() )
        if not omega(mol):
            print( "omega failed on %s" % mol.GetTitle() )
        title= mol.GetTitle()
        oechem.OETriposAtomNames( mol)
        oequacpac.OEAssignPartialCharges(mol, oequacpac.OECharges_AM1BCCSym)
        oechem.OEWriteConstMolecule(ofs, mol)
    
    ofs.close()
    return 0


#############################################################################
# INTERFACE
#############################################################################

InterfaceData = '''
!CATEGORY "input/output options"
    !PARAMETER -in
      !ALIAS -i
      !TYPE string
      !REQUIRED true
      !KEYLESS 1
      !VISIBILITY simple
      !BRIEF Input filename
    !END

    !PARAMETER -out
      !ALIAS -o
      !TYPE string
      !REQUIRED true
      !KEYLESS 2
      !VISIBILITY simple
      !BRIEF Output filename
    !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
