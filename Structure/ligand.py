#!/usr/bin/env python
'''
This module defines the class Ligand which implements
attributes and methods required to handle ligand (small molecules).
'''
import sys
import numpy
from molecule import Molecule, Atom

class Ligand(Molecule, Atom):
    'Base class for a Ligand'
    
    def __init__(self):
        Molecule.__init__(self)
        self.molid = Molecule.molid
        self.nor = 0   # Number of residues in a ligand; WILL ALWAYS BE 1
        self.resids = [] # List of residue ids of a ligand
        self.atmidx = [] # List of atom indices from the input file 
        self.residue = {}   # Information about each residue. Key: Residue id.
    
    def molecule_type(self):
        'Return a string representing the type of molecule this is.'
        return "Ligand"