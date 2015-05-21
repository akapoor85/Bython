#!/usr/bin/env python

from abc import ABCMeta, abstractmethod

class Molecule:
    'Abstract Base class for a molecule'
    
    __metaclass__ = ABCMeta
    molid = 0 # Bython molid to keep internal count of number of molecules
        
    def __init__(self):
        Molecule.molid +=1
        
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
    

    

        
    