#!/usr/bin/env python
'''
This module defines central functions such as 'Loadstruc' to load
molecule structure file. Loadstruc further makes appropriate function calls
depending on the type of input file to be loaded.

Molecule can also be directly loaded from www.pdb.org by providing PDB id
of the molecule to function FetchPDB.
'''
import sys
import os
import urllib
from configstruc import *
from loadpdbmulti import Loadpdb

def FetchPDB(pdbid=None,destination=None):
    '''
    Download pdb file corresponding to pdbid from www.pdb.org.
    The file is saved in the folder specified by destination with
    the name pdbid.pdb
    '''
    try:
        assert(pdbid != None) #Check if pdbid was passed
    except AssertionError:
        sys.exit("Provide a PDB id to download")
    
    url= 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
       
    if destination == None:
        destination = os.getcwd()
        
    filename= destination+"/"+pdbid+".pdb"
    urllib.urlretrieve(url,filename)

def Loadstruc(intype=None, moltype=None, filename=None, hetatm = 'yes'):
    '''
    Load an input molecule structure file.
    intype: filetype of input file (ex: PDB, mol2 etc.)
    moltype: Molecule type (ex: Protein, Lipid, ligand etc.)
    hetatm: Load data for hetatm records (default: yes)
    '''
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
        return Loadpdb(pdb=strucfile, molpdb=moltype, hetatm=hetatm)
    

        
    
    