#!/usr/bin/env python
'''
Each HETATM residue is set as a separate Ligand object
'''
import sys
import numpy
from copy import deepcopy
from molecule import Molecule, Atom
from protein import Protein
from ligand import Ligand
from lipid import Lipid
from configstruc import Mol_types

def Pdbcordsec(atomline):
    'Given a atom coordinate section line from PDB file return information about atom'
    AtomName= atomline[12:16].upper().strip()
    ResName= atomline[17:21].upper().strip() #PDB only allows 3-letter resname but here upto four letters will work
    Chain= atomline[21].upper().strip()
    ResNo= int(atomline[22:26].upper().strip())
    CordX= float(atomline[30:38].upper().strip())
    CordY= float(atomline[38:46].upper().strip())
    CordZ= float(atomline[46:54].upper().strip())
    Occ  = float(atomline[54:60].upper().strip())
    Bfac = float(atomline[60:66].upper().strip())
    return AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac   

def Loadpdb(pdb=None, hetatm= True, verbose=False):
    try:
        assert(pdb != None) #Check if filehandle to PDB file is passed
    except AssertionError:
        sys.exit("**No filehandle passed**. Pass a filehandle (to a pdb file) as an argument to Loadpdb.")
    
    AtomNumber=0 #Keeps track of atom indices (assigned in the order atoms listed in input file)
    mol_data={} #Key: molid; Value: Molecule_Type Object; Keep track of different molecules (different chains or molecule type) in input structure
    first_res =True  #To identify molecule type of every molecule in input structure and accordingly define Molecule object.
    Prev_res=0 # to keep track of residue change 
    Prev_chain='aa' # to keep track of chain change in HETATM section
    atmTohet = True #To determine transition from ATOM to HETATM record
    frame_tag = '' # To keep track of multi-frame entry (multiple entry for same molecule type with same chain id)
    '''Load the PDB structure file'''
    for line in pdb:
        if line[0:4]=="ATOM" and line[12:16].upper().strip() not in ["OXT"]:
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            
            #Check for unrecognized residue and new molecule
            if not first_res:
                if ResName.lower() not in Mol_types['protein'] + Mol_types['lipid'] + Mol_types['ligand']:
                    print "*** Unrecognized residue name: "+ ResName+ " ***."
                    sys.exit("In file configstruc.py: Add missing residue name("+ ResName+ ") to appropriate molecule in Mol_types")
                elif ResNo != Prev_res or (Prev_chain != Chain and mol.molecule_type().lower()=='ligand'):
                    if mol.molecule_type().lower()=='ligand':
                        mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
                        first_res = True
                    elif Prev_chain != Chain or ResName.lower() not in Mol_types[mol.molecule_type().lower()]: #Either Chain is different or New residue doesn't belong to current molecule type
                        mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
                        first_res = True
                    elif frame_tag.lower() in ['endmdl', 'ter', 'end'] and Prev_chain == Chain: #Different molecule (of same molecule type) with same chain id; as in trajectory frames
                        mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
                        first_res = True
                        
                        
            if first_res: #Initialize mol for new chain or molecule
                if ResName.lower() in Mol_types['protein']:
                    mol = Protein()
                elif ResName.lower() in Mol_types['lipid']:
                    mol = Lipid()
                elif ResName.lower() in Mol_types['ligand']:
                    mol = Ligand()
                else:
                    print "*** Unrecognized residue name: "+ ResName+ " ***.\n Cannot initialize Molecule object."
                    sys.exit("In file configstruc.py: Add missing residue name("+ ResName+ ") to appropriate molecule in Mol_types")
                first_res = False
                frame_tag = ''
            
            mol.AddToResidue(atom=atm, occ=Occ, bfac=Bfac)
            mol.atmidx.append(AtomNumber)
            Prev_res= ResNo
            Prev_chain=Chain
        elif line[0:6]=="HETATM" and hetatm == True:
            if atmTohet:
                mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
                first_res = True
                atmTohet = False
            AtomNumber+=1
            AtomName, ResName, Chain, ResNo, CordX, CordY, CordZ, Occ, Bfac = Pdbcordsec(line)
            atm = Atom(AtomName, AtomNumber, ResName, Chain, ResNo, CordX, CordY, CordZ)
            
            #Check for new ligand molecule 
            if not first_res and (ResNo != Prev_res or Prev_chain != Chain):
                mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
                first_res = True
            
            #Initialize mol for new chain or molecule            
            if first_res: 
                if ResName.lower() in Mol_types['ligand']:
                    mol = Ligand()
                else:
                    print "*** Unrecognized residue name: "+ ResName+ " ***.\n Cannot initialize Ligand object."
                    sys.exit("In file configstruc.py: Add missing residue name ("+ ResName+ ") to ligand molecule in Mol_types")
                first_res = False
            
            mol.AddToResidue(atom=atm, occ=Occ, bfac=Bfac)
            mol.atmidx.append(AtomNumber)
            Prev_res= ResNo
            Prev_chain=Chain
        elif line[0:3].lower() in ["ter", "end"] or line[0:6].lower() == "endmdl":
            frame_tag = line[0:3]
    #append the last mol object to mol_data
    mol_data[Molecule.molid] = deepcopy(mol) #copy mol object into dictionary
    
    if verbose:
        print "Number of molecules in input file: ", len(mol_data), "\n"
    #Update mol_data[molid].nor, mol_data[molid].resids, and check for chain breaks in non-ligand molecules
    for key in sorted(mol_data):
        if verbose:
            print "Molid:", key,"Molecule_Type:",mol_data[key].molecule_type()
        mol_data[key].resids = sorted(mol_data[key].residue)
        mol_data[key].nor = len(mol_data[key].resids)
        if mol_data[key].molecule_type().lower() != 'ligand':
            #Check for chain breaks in non-ligand molecule
            resids_diff=numpy.array(mol_data[key].resids[1:]) - numpy.array(mol_data[key].resids[:-1])
            if mol_data[key].nor != (numpy.sum(resids_diff)+1):
                break_indices = numpy.where(resids_diff > 1)
                print "Chain break encountered in molecule",key, "at residue positions: "
                for res in break_indices[0]:
                    print mol_data[key].resids[res],
                print "\n"
                mol_data[key].chain_break = True
    return mol_data 