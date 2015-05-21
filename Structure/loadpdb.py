#!/usr/bin/env python

import sys
import numpy
from molecule import Atom
from protein import Protein
from configstruc import Mol_types

def Pdbcordsec(atomline):
    'Given a atom coordinate section line from PDB file return information about atom'
    AtomName= atomline[12:16].upper().strip()
    ResName= atomline[17:20].upper().strip()
    Chain= atomline[21].upper().strip()
    ResNo= int(atomline[22:26].upper().strip())
    CordX= float(atomline[30:38].upper().strip())
    CordY= float(atomline[38:46].upper().strip())
    CordZ= float(atomline[46:54].upper().strip())
    Occ  = float(atomline[54:60].upper().strip())
    Bfac = float(atomline[60:66].upper().strip())
    return AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac   

def Loadpdb(pdb=None, molpdb='protein', hetatm= 'yes'):
    try:
        assert(pdb != None) #Check if filehandle to PDB file is passed
    except AssertionError:
        sys.exit("**No filehandle passed**. Pass a filehandle (to a pdb file) as an argument to Loadpdb. ")
    
    try:
        assert(molpdb.lower() in Mol_types) #Check if molecule type recognized
        '''Load a PDB structure file of a molecule of type molpdb'''
        if molpdb.lower() == 'protein':
            mol = Protein()
    except AssertionError:
        sys.exit("Bython does not handle "+molpdb+" molecule")
    
    AtomNumber=0 # For Assigning Atom Index
    
    for line in pdb:
        if line[0:4]=="ATOM" and line[12:16].upper().strip() not in ["OXT"]:
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            mol.AddToResidue(atom=atm, occ=Occ, bfac=Bfac, seg='Protein')
            mol.resids.append(ResNo)
        elif line[0:6]=="HETATM" and hetatm.lower() == 'yes':
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            mol.AddToResidue(atom=atm, occ=Occ, bfac=Bfac, seg='Hetero')
    #Remove duplicate entries from mol.resid, and update mol.nor.
    #*Might need to put under if moltype protein in future*
    mol.resids = list(set(mol.resids))
    mol.nor = len(mol.resids)
    #Check for chain breaks in protein
    resids_diff=numpy.array(mol.resids[1:]) - numpy.array(mol.resids[:-1])
    if mol.nor != (numpy.sum(resids_diff)+1):
        break_indices = numpy.where(resids_diff > 1)
        print "Chain breaks encountered at residue positions: "
        for res in break_indices[0]:
            print mol.resids[res],
        
    return mol