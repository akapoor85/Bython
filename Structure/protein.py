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
from configstruc import RESIDUE_CHI_CHARMM

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
    
    def GetChi(self, resid=None, atom_names= 'charmm'):
        '''
        Returns all chi dihedrals for a given residue id. Atom names specifies the force-field for atom names in the input file.
        Return type: Dict; Key: Chi1...Chi5 (depending on residue); Value: [[chi_atom_indices], chi_dihedral_value_in_degrees]
        '''
        VALID_FF = ['charmm']
        try:
            assert(resid != None) #Check if resid is passed
        except AssertionError:
            raise ValueError('*** No residue id provided ***')
        
        try:
            assert(atom_names.lower() in VALID_FF) 
        except AssertionError:
            raise ValueError('*** Unrecognized atom_names parameter value %s *** \nValid inputs for atom_names: %s' % (atom_names, ','.join(VALID_FF)))
        
        #Built atom_name -> atom_index map for the given residue
        atom_map = {self.residue[resid]['atoms'][atmidx]['name'] : atmidx for atmidx in self.residue[resid]['atoms']}
        
        if atom_names.lower() == 'charmm':
            chi_dict = { 'chi'+str(chi_no+1) : [[atom_map[dihed[0]], atom_map[dihed[1]], atom_map[dihed[2]], atom_map[dihed[3]]
                                                 ], self.GetDihed(atom_map[dihed[0]], atom_map[dihed[1]], atom_map[dihed[2]], atom_map[dihed[3]])]
                for chi_no, dihed in enumerate(RESIDUE_CHI_CHARMM[self.residue[resid]['name']])}
                   
        return chi_dict        
    
    def molecule_type(self):
        'Return a string representing the type of molecule this is.'
        return "Protein"