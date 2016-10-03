#!/usr/bin/env python
import os

def Chi_dict(chi_file):
    '''
    Reads a file with information about chi dihedrals for each residue and returns the corresponding dict.
    '''
    res_chi={}
    for line in chi_file:
        if line.strip(): #skip blank lines
            line = line.strip().split('\t')
            if not res_chi.has_key(line[0].upper()):
                res_chi[line[0].upper()] = [entry.upper().split('-') for entry in line[1:]]
    return res_chi            
         
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

#Dictionary of Chi dihedrals for every protein residue
file_res_chi_charmm = open(os.path.join(os.path.dirname(__file__),'residue_chi_charmm'), 'r')
headerline = file_res_chi_charmm.readline()
RESIDUE_CHI_CHARMM = Chi_dict(file_res_chi_charmm)


        

