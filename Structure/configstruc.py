#!/usr/bin/env python

Struc_in_formats = ['pdb'] #Recognized Structure File Types; MUST be in lowercase
#Recognized Molecule Types; MUST be in lowercase
Mol_types = {'protein': ['ala', 'arg', 'asn', 'asp', 'cys', 'glu',
                         'gln', 'gly', 'his', 'ile', 'leu', 'lys',
                         'met', 'phe', 'pro', 'ser', 'thr', 'trp',
                         'tyr', 'val', 'hid', 'hie', 'hsd', 'ace',
                         'nma'],
            'lipid': ['pop', 'popc'],
            'ligand': ['suv', 'ola', 'hoh', 'hem', 'chl', 'chl1', 'nal', 'jdc',
                       'cit', 'olc', 'peg', 'mg', 'gnp']} 