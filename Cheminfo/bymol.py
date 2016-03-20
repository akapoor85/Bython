#!/usr/bin/env python

from __future__ import division
from ChemUtil import loadmol
import sys


class Bymol:
    '''
    Base class for small molecules used in cheminformatic analysis.
    Since RDKit is used, Bymol has access to all RDKit mol attributes.
    '''
    count = 0 # Keeps track of number of molecules
    def __init__(self):
        self.by_mol = None #holds rdkit mol object
        self.molid = Bymol.count + 1
        self.name = ''
        self.atoms = None #holds RDKit Atom Sequence object
        self.fpt = None #Holds a user-requested 2D fingerprint
        self.property = {} #To set user-defined properties for a molecule
        Bymol.count = Bymol.count + 1
    
    def Readmol(self, intype=None, molstring=None, name='',remH = False):
        '''
        Reads a small molecule file (molstring) of type (intype: pdb, mol2, mol) or a SMILES
        string (intype: smi) and assigns the returned RDKit mol object
        to Bymol instance attribute by_mol. A molecule name can be assigned wih name. 
        By default, RDKit removes hydrogens (remH) when reading pdb or mol2 file. Here, we will
        retain all explicit hydrogens in structure files. In SMILES, hydrogens are treated implicit even
        when defined explicitly. 
        '''
        valid_intype = ['mol2', 'pdb', 'smi', 'mol']
        try:
            assert(intype.lower() in valid_intype) #Check if input format recognized
        except AssertionError:
            sys.exit("Invalid intype argument: "+intype+".\nChoose one of %s" % (','.join(valid_intype)))
        
        #Load the mol using RDKit functionality
        self.by_mol = loadmol(intype=intype, molstring=molstring, remH=remH)
        self.name = name
        self.atoms = self.by_mol.GetAtoms()
        