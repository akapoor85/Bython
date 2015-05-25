#!/usr/bin/env python
'''
ATOM record must always end in TER even if there is just one macromolecule (i.e. no multiple chains)
ATOM record of multiple chains must be separated by TER
ATOM and HETATM must be separated by TER
'''
import sys
import numpy
from copy import deepcopy
from molecule import Molecule, Atom
from protein import Protein
from ligand import Ligand
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
        mol_data={} #Key: molid; Keep track of different molecules (different chains) in a input structure
        '''Load a PDB structure file of a molecule of type molpdb'''
        if molpdb.lower() == 'protein':
            mol = Protein()
    except AssertionError:
        sys.exit("Bython does not handle "+molpdb+" molecule")
    
    AtomNumber=0 # For Assigning Atom Index
    check_mol =False #To keep track of new chain
    check_het =False #To keep track of new Hetero residue
    Prev_res=0 # to keep track of residue change in HETATM section; a new Molecule object is assigned for every residue.
    Prev_chain='a' # to keep track of chain change in HETATM section
    for line in pdb:
        if check_mol:
            mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
            if line[0:4]=="ATOM":
                if molpdb.lower() == 'protein': #reinitialize mol for new chain
                    mol = Protein()
            elif line[0:6]=="HETATM" and hetatm.lower() == 'yes':
                mol = Ligand()
            check_mol = False
        if line[0:4]=="ATOM" and line[12:16].upper().strip() not in ["OXT"]:
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            mol.AddToResidue(atom=atm, occ=Occ, bfac=Bfac, seg='Protein')
            mol.resids.append(ResNo)
        elif line[0:3]=="TER":
            check_mol = True #A new chain will be added if another chain with ATOM record follows
        elif line[0:6]=="HETATM" and hetatm.lower() == 'yes':
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            if (ResNo != Prev_res or Prev_chain != Chain) and check_het == True: #For first HETATM check_het is always False
                mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
                mol = Ligand()
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            mol.AddToResidue(atom=atm, occ=Occ, bfac=Bfac, seg='Hetero')
            if Prev_res == 0:
                check_het = True 
            Prev_res= ResNo
            Prev_chain=Chain
    if hetatm.lower() == 'yes': #HETATM record was added; append the last hetero residue object to mol_data
        mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
    
    print "Number of molecules in input file: ", len(mol_data)
    #Remove duplicate entries from mol_data[molid].resids, and update mol_data[molid].nor
    for key in sorted(mol_data):
        print "\nmolid:", key, mol_data[key].molecule_type()
        if mol_data[key].molecule_type().lower() != 'ligand':
            mol_data[key].resids = list(set(mol_data[key].resids))
            mol_data[key].nor = len(mol_data[key].resids)
            #Check for chain breaks in protein
            resids_diff=numpy.array(mol_data[key].resids[1:]) - numpy.array(mol_data[key].resids[:-1])
            if mol_data[key].nor != (numpy.sum(resids_diff)+1):
                break_indices = numpy.where(resids_diff > 1)
                print "Chain break encountered in molecule",key, "at residue positions: "
                for res in break_indices[0]:
                    print mol_data[key].resids[res],
        
    return mol_data 