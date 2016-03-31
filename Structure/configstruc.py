#!/usr/bin/env python

Struc_in_formats = ('pdb') #Recognized Structure File Types by Structure sub-package; MUST be in lowercase
#Recognized Molecule Types; MUST be in lowercase
Mol_types = {'protein': ('ala', 'arg', 'asn', 'asp', 'cys', 'glu',
                         'gln', 'gly', 'his', 'ile', 'leu', 'lys',
                         'met', 'phe', 'pro', 'ser', 'thr', 'trp',
                         'tyr', 'val', 'hid', 'hie', 'hsd', 'ace',
                         'nma'),
            'lipid': ('pop', 'popc', 'chl', 'chl1'),
            'ligand': ('suv', 'ola', 'hoh', 'hem', 'nal', 'jdc',
                       'cit', 'olc', 'peg', 'mg', 'gnp', 'gdp',
                       'gtp', 'atp', 'adp', 'mtd', 'ms1'),
            'ion': ('mg', 'cl', 'na', 'ca', 'zn')}

Mol_types['ligand'] = Mol_types['ligand'] + Mol_types['ion']

#List of atoms symbol considered to be non-polar and polar
NON_POLAR_ATOMS = ('c')
POLAR_ATOMS = ('n', 'o', 's')

#List of protein aromatic residues
PROTEIN_AROMATIC_RES = ('his', 'hid', 'hie', 'hsd', 'trp', 'tyr', 'phe')
#List of protein charged residues
PROTEIN_POSITIVE = ('arg', 'lys', 'hid', 'hie', 'hsd', 'his')
PROTEIN_NEGATIVE = ('asp', 'glu')