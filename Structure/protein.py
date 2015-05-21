#!/usr/bin/env python
'''
This module defines the class Protein which implements
attributes and methods required to handle protein macromolecule.

The basic building block is a single residue. All the information
about each residue in a protein is stored in a nested dictionary "residue".
The structure of dictionary residue is:
residue[key-> residue id] = {name: resname,
                             chain: chain id,
                             segment: segment name (->either Protein or Hetero),
                             atoms:{atom number: {name: atom name, cord: [cordx, cordy, cordz],
                                                  B: b-factor, O: occupancy}
                                    }
                            }
Each residue dictionary has a nested dictionary 'atoms' that contains data
related to each atom in a residue. The key to dictionary atoms is atom number
(numbered in the order they appear in input file).
'''
import sys
import numpy
from molecule import Molecule, Atom

class Protein(Molecule, Atom):
    'Base class for a Protein'
    
    def __init__(self):
        Molecule.__init__(self)
        self.nor = 0   # Number of residues in a protein
        self.nohet = 0 # Number of hetero residues in a protein
        self.chains = [] # List of different chains of a protein
        self.resids = [] # List of residue ids of a protein
        self.residue = {}   # Information about each residue. Key: Residue id.
    
    def AddToResidue(self, atom=None, seg=None, occ=1.0, bfac=0.0):
        '''
        Add atom information to build a residue (one atom at a time).
        Argument atom is the class Atom object.
        occ and bfac correspond to occupancy and b-factor of an atom
        represented by object atom.
        seg is the segment (Protein or Hetero) to which an atom belongs.
        '''
        try:
            assert(isinstance(atom, Atom)) #Check if atom is an object of class Atom
        except AssertionError:
            sys.exit('** Argument atom must be object of class Atom **.\nInvalid atom argument to AddToResidue.')
        
        try:
            assert(seg.lower() in ['protein', 'hetero']) #Check if seg is passed correctly
        except AssertionError:
            sys.exit('** Argument seg must be either Protein or Hetero **.\nInvalid seg argument to AddToResidue.')
        
        if not self.residue.has_key(atom.resno):
            self.residue[atom.resno]= {'name': atom.resname, 'chain': atom.chain, 'segment': seg,
                                       'atoms':{atom.idx:{'name':atom.name, 'cord': [atom.cordx,atom.cordy,atom.cordz],
                                        'B': bfac, 'O': occ}}}
        else:
            self.residue[atom.resno]['atoms'][atom.idx]= {'name':atom.name,'cord': [atom.cordx,atom.cordy,atom.cordz],'B': bfac, 'O': occ}
        
    def GetNor(self):
        'Returns number of residues in the segment protein'
        self.nor = len([key for key in self.residue if self.residue[key]['segment'].lower()=='protein'])
        return self.nor
    
    def GetNohet(self):
        'Returns number of hetero residues in a protein (from segment Hetero)'
        self.nohet = len([key for key in self.residue if self.residue[key]['segment'].lower()=='hetero'])
        return self.nohet
    
    def GetChains(self):
        'Returns id of different chains in a protein'
        self.chains = list(set([self.residue[key]['chain'] for key in self.residue if self.residue[key]['segment'].lower()=='protein']))
        return self.chains
    
    def GetCord(self, resid= None, hetatm = 'no'):
        '''
        Returns coordinates (as numpy array) of all atoms (ordered by atom index) of a
        residue specfied by residue id or for all atoms in a protein (resid: all)
        By default only segment protein (hetatm = 'no') is returned for resid 'all'. 
        '''
        try:
            assert(resid != None) #Check if resid is passed
        except AssertionError:
            sys.exit('** No residue id provided **. \nValid input to resid: residue id or all.')
            
        if str(resid).lower() == 'all' and hetatm.lower() == 'yes':
            resid = sorted(self.residue)
        if str(resid).lower() == 'all' and hetatm.lower() == 'no': # Extract id for segment Protein only
            resid = [key for key in sorted(self.residue) if self.residue[key]['segment'].lower() == 'protein']
        else:
            resid = [resid]
                  
        coordinates = [self.residue[key]['atoms'][atmidx]['cord'] for key in resid for atmidx in sorted(self.residue[key]['atoms'])]
        return numpy.array(coordinates)
    
    def GetResids(self):
        '''
        Returns a list of residue ids in a protein. Hetero segment is not included.
        For residue ids of protein+hetero, just use Objectname.residue.keys()
        '''
        self.resids = [key for key in sorted(self.residue) if self.residue[key]['segment'].lower() == 'protein']
        return self.resids
        
    def molecule_type(self):
        'Return a string representing the type of molecule this is.'
        return "Protein"