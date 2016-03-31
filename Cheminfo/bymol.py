#!/usr/bin/env python

from __future__ import division
from ChemUtil import Loadmol
import sys
import fingerprints


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
        self.property = {} #To set default and user-defined properties for a molecule. key: atom index, value: {property name: property value}
        self.neighbors = {} #For each atom, a tuple of directly bonded neighbors. key: atom index, value: (list of neighbor indices)
        self.rings = () #A tuple of all rings (atom-indices)
        self.n_aromatic = 0
        self.n_heavy = 0
        Bymol.count = Bymol.count + 1
    
    def Readmol(self, intype=None, molstring=None, name='',remH = False, no_sani=False):
        '''
        Reads a small molecule file (molstring) of type (intype: pdb, mol2, mol) or a SMILES
        string (intype: smi) and assigns the returned RDKit mol object
        to Bymol instance attribute by_mol. A molecule name can be assigned wih name. 
        By default, RDKit removes hydrogens (remH) when reading pdb or mol2 file. Here, we will
        retain all explicit hydrogens in structure files. In SMILES, hydrogen treatment is implicit even
        when defined explicitly. Set no_sani to False (no sanitization) to turn off RDKit sanitization. 
        '''
        valid_intype = ['mol2', 'pdb', 'smi', 'mol']
        try:
            assert(intype.lower() in valid_intype) #Check if input format recognized
        except AssertionError:
            sys.exit("Invalid intype argument: "+intype+".\nChoose one of %s" % (','.join(valid_intype)))
        
        #Load the mol using RDKit functionality
        self.by_mol = Loadmol(intype=intype, molstring=molstring, remH=remH, no_sani=no_sani)
        #Set different properties
        self.name = name
        self.atoms = self.by_mol.GetAtoms()
        #Set neighbors and property dicts
        count_aro = 0 #Number of aromatic atoms
        for atom in self.by_mol.GetAtoms():
            for neighb in atom.GetNeighbors():
                if not self.neighbors.has_key(atom.GetIdx()):
                    self.neighbors[atom.GetIdx()] = (neighb.GetIdx(),)
                else:
                    self.neighbors[atom.GetIdx()]+= (neighb.GetIdx(),)
            self.property[atom.GetIdx()] = {'symbol': atom.GetSymbol()}
            self.property[atom.GetIdx()]['aromatic'] = atom.GetIsAromatic()
            if self.property[atom.GetIdx()]['aromatic']:
                count_aro +=1
            self.property[atom.GetIdx()]['inring'] = atom.IsInRing()
            self.property[atom.GetIdx()]['charge'] = 0.0
        ring=self.by_mol.GetRingInfo()
        self.rings = ring.AtomRings()
        self.n_aromatic = count_aro
        self.n_heavy = self.by_mol.GetNumHeavyAtoms()
    
    def Set_Fingerprint(self, fpname='linear', bitsize=2048, fptype='bit', diameter=4):
        '''
        Assigns the 2D fingerprint to class instance's fpt attribute.
        INPUT-> fpname: Name of fingerprint to calculate (one of linear, maccs,
        atompairs, torsion, ecfp, and fcfp)
        bitsize: bit vector length (not used when fpname maccs)
        fptype: one of bit or count (only valid for a subset of fpname)
        '''
        valid_fpname = ('linear', 'maccs', 'atompairs', 'torsion', 'ecfp', 'fcfp')
        valid_fpname_count = ('atompairs', 'torsion','ecfp', 'fcfp')
        #Validate Input
        try:
            assert(fpname.lower() in valid_fpname)
        except AssertionError:
            sys.exit("***ERROR: Unrecognized fingerprint type (%s) requested.***\n\
                     Valid options are %s" % (fpname, ', '.join(valid_fpname)))
            
        if fptype.lower() == 'count' and fpname.lower() not in valid_fpname_count:
            print "%s fingerprint does not return %s, reverting to fptype=bit" % (fpname.upper(), fptype)
            fptype='bit'
        elif fptype.lower() not in ('bit', 'count'):
            print "Unrecognized fptype passed: %s, reverting to fptype=bit" % fptype
            fptype='bit'
        
        if fpname.lower() == 'torsion' and fptype.lower() == 'bit':
            print "%s fingerprint does not return %s, reverting to fptype=count" % (fpname.upper(), fptype)
            fptype='count'
        
        #Call appropriate function
        if fpname.lower() == 'linear':
            self.fpt = fingerprints.Fptlinear(self.by_mol, fpSize=bitsize)
        elif fpname.lower() == 'maccs':
            self.fpt = fingerprints.FptMACCS(self.by_mol)
        elif fpname.lower() == 'atompairs':
            self.fpt = fingerprints.FptAtomPairs(self.by_mol, fptype=fptype)
        elif fpname.lower() == 'torsion':
            self.fpt = fingerprints.FptTorsion(self.by_mol)
        elif fpname.lower() == 'ecfp':
            self.fpt = fingerprints.FptCircular(self.by_mol, radius=diameter/2, nBits=bitsize, fptype=fptype)
        elif fpname.lower() == 'fcfp':
            self.fpt = fingerprints.FptCircular(self.by_mol, radius=diameter/2, nBits=bitsize, fptype=fptype, feature=True)
            
            
        
        
        
        
        