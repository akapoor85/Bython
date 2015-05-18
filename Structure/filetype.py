#!/usr/bin/env python

from molecule import *
from protein import Protein

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

def Loadpdb(pdb, molpdb):
    '''Load a PDB structure file of a molecule of type molpdb'''
    if molpdb.lower() == 'protein':
        mol = Protein()
    
    AtomNumber=0 # For Assigning Atom Index
    for line in pdb:
        if line[0:4]=="ATOM" and line[12:16].upper().strip() not in ["OXT"]:
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            mol.AddResidue(atm, Occ, Bfac)
            
    return mol