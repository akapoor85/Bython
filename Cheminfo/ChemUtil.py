#!/usr/bin/env python

from __future__ import division
from rdkit import Chem
import sys
import sanifix4

def loadmol(intype=None, molstring=None, remH = False):
    #First load input mol as RDKit mol without sanitizing
    if intype.lower() == 'pdb':
        temp_mol = Chem.MolFromPDBFile(molstring, removeHs=remH, sanitize=False)
    elif intype.lower() == 'mol2':
        temp_mol = Chem.MolFromMol2File(molstring, removeHs=remH, sanitize=False)
    elif intype.lower() == 'mol':
        temp_mol = Chem.MolFromMolFile(molstring, removeHs=remH, sanitize=False)
    elif intype.lower() == 'smi':
        temp_mol = Chem.MolFromSmiles(molstring, sanitize=False)
    else:
        sys.exit("Invalid intype argument: "+intype+".\nChoose one of %s" % (','.join(valid_intype)))
        
    #Sanitize mol
    try:
        Chem.SanitizeMol(temp_mol)
        #Sanitization successful when no exception raised
        return temp_mol
    except ValueError as err:
        if "can't kekulize" in err.message.lower():
            #Attempt to sanitize it with sanifix4; need to reload molecule as SanitizeMol does some modifications
            if intype.lower() == 'pdb':
                temp_mol_new = Chem.MolFromPDBFile(molstring, removeHs=remH, sanitize=False)
            elif intype.lower() == 'mol2':
                temp_mol_new = Chem.MolFromMol2File(molstring, removeHs=remH, sanitize=False)
            elif intype.lower() == 'mol':
                temp_mol_new = Chem.MolFromMolFile(molstring, removeHs=remH, sanitize=False)
            elif intype.lower() == 'smi':
                temp_mol_new = Chem.MolFromSmiles(molstring, sanitize=False)
            
            temp_mol_fixed = sanifix4.AdjustAromaticNs(temp_mol_new)
            if temp_mol_fixed is not None: #sanitization was successful
                return temp_mol_fixed
            else:
                sys.exit('***ERROR: Cannot kekulize molecule %s. SANIFIX4 did not work either.***' % molstring)
        else:
            sys.exit('***ERROR: RDKit failed sanitization for molecule %s***\n%s' % (molstring, err.message))
            
            
    
    
        
    
