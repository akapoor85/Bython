#!/usr/bin/env python
'''
This module defines the class Protein which implements
attributes and methods required to handle protein macromolecule.

The basic building block is a single residue. All the information
about each residue in a protein is stored in a nested dictionary "residue".
The structure of dictionary residue is:
residue[key-> residue id] = {name: resname,
                             chain: chain id,
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
        self.molid = Molecule.molid
        self.nor = 0   # Number of residues in a protein
        self.resids = [] # List of residue ids of a protein
        self.atmidx = [] # List of atom indices  
        self.residue = {}   # Information about each residue. Key: Residue id.
        self.chain_break = False # This will be True if chain breaks are encountered
    
    def molecule_type(self):
        'Return a string representing the type of molecule this is.'
        return "Protein"