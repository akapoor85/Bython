#!/usr/bin/env python

from molecule import *

class Protein(Molecule, Atom):
    'Base class for a Protein'
    
    def __init__(self):
        Molecule.__init__(self)
        self.nor = 0   # Number of residues
        self.residue = {}   # Information about each residue. Key: Resid
    
    def AddResidue(self, atom, occ, bfac):
        '''Add atoms to build residue of a protein. atom is the class Atom object.
        Occ and Bfac correspond to occupancy and b-factor of atom
        represented by object atom'''
        
        if not self.residue.has_key(atom.resno):
            self.residue[atom.resno]= {'name': atom.resname, 'chain': atom.chain,
                                       'atoms':{atom.idx:{'name':atom.name, 'cord': [atom.cordx,atom.cordy,atom.cordz],
                                        'B': bfac, 'O': occ}}}
        else:
            self.residue[atom.resno]['atoms'][atom.idx]= {'name':atom.name,'cord': [atom.cordx,atom.cordy,atom.cordz],'B': bfac, 'O': occ}
        
    
    def molecule_type(self):
        """"Return a string representing the type of molecule this is."""
        return "Protein"