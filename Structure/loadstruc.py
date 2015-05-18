#!/usr/bin/env python

from configstruc import *
import filetype as FileT
import sys

def Loadstruc(intype=None, moltype=None, filename=None):
    '''Read a molecule structure file of type intype containing a molecule of type moltype'''
    try:
        assert(moltype.lower() in Mol_types) #Check if molecule type recognized
    except AssertionError:
        sys.exit("Bython does not handle "+moltype+" molecule")
    
    try:
        assert(intype.lower() in Struc_in_formats) #Check if input format recognized
        strucfile = open(filename, 'r')
    except AssertionError:
        sys.exit("Bython does not handle "+intype+" input files")
    except IOError:
        sys.exit("Cannot locate file "+ filename+ ". Invalid filename or path!")
    
    if intype.lower() == 'pdb':
        return FileT.Loadpdb(strucfile, moltype)
    

        
    
    