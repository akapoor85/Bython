#!/usr/bin/env python

from abc import ABCMeta, abstractmethod

class Molecule:
    'Abstract Base class for a molecule'
    
    __metaclass__ = ABCMeta
    molid = 0 # Bython molid to keep internal count of number of molecules
        
    def __init__(self):
        Molecule.molid +=1
        
    def AddToResidue(self, atom=None, seg=None, occ=1.0, bfac=0.0):
        '''
        Add atom information to build a residue (one atom at a time).
        Argument atom is the class Atom object.
        occ and bfac correspond to occupancy and b-factor of an atom
        represented by object atom.
        seg is the segment (Protein or Hetero) to which an atom belongs.
        '''
        try:
            assert(isinstance(atom, Atom)) #Check if atom is an object of class Atom
        except AssertionError:
            sys.exit('** Argument atom must be object of class Atom **.\nInvalid atom argument to AddToResidue.')
        
        try:
            assert(seg.lower() in ['protein', 'hetero']) #Check if seg is passed correctly
        except AssertionError:
            sys.exit('** Argument seg must be either Protein or Hetero **.\nInvalid seg argument to AddToResidue.')
        
        if not self.residue.has_key(atom.resno):
            self.residue[atom.resno]= {'name': atom.resname, 'chain': atom.chain, 'segment': seg,
                                       'atoms':{atom.idx:{'name':atom.name, 'cord': [atom.cordx,atom.cordy,atom.cordz],
                                        'B': bfac, 'O': occ}}}
        else:
            self.residue[atom.resno]['atoms'][atom.idx]= {'name':atom.name,'cord': [atom.cordx,atom.cordy,atom.cordz],'B': bfac, 'O': occ}
    
    @abstractmethod
    def molecule_type(self):
        """"Return a string representing the type of molecule this is."""
        pass
    
class Atom:
    'Base class for an atom'
    
    def __init__(self, name, idx, resN, chain, resNo, cordx, cordy, cordz):
        self.name = name
        self.idx  = idx
        self.resname  = resN
        self.chain  = chain
        self.resno  = resNo
        self.cordx  = cordx
        self.cordy  = cordy
        self.cordz  = cordz
    

    

        
    