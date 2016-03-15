#!/usr/bin/env python

from __future__ import division
from abc import ABCMeta, abstractmethod
from copy import deepcopy
import numpy
import sys
import os

class Molecule:
    'Abstract Base class for a molecule'
    
    __metaclass__ = ABCMeta
    molid = 0 # Bython molid to keep internal count of number of molecules
        
    def __init__(self):
        Molecule.molid +=1
        
    def AddToResidue(self, atom=None, occ=1.0, bfac=0.0):
        '''
        Add atom information to build a residue (one atom at a time).
        Argument atom is the class Atom object.
        occ and bfac correspond to occupancy and b-factor of an atom
        represented by object atom.
        '''
        try:
            assert(isinstance(atom, Atom)) #Check if atom is an object of class Atom
        except AssertionError:
            sys.exit('** Argument atom must be object of class Atom **.\nInvalid atom argument to AddToResidue.')
        
        if not self.residue.has_key(atom.resno):
            self.residue[atom.resno]= {'name': atom.resname, 'chain': atom.chain,
                                       'atoms':{atom.idx:{'name':atom.name, 'cord': [atom.cordx,atom.cordy,atom.cordz],
                                        'B': bfac, 'O': occ}}}
        else:
            self.residue[atom.resno]['atoms'][atom.idx]= {'name':atom.name,'cord': [atom.cordx,atom.cordy,atom.cordz],'B': bfac, 'O': occ}
    
    def GetCord(self, resid= None):
        '''
        Returns coordinates (as numpy array) of all atoms (ordered by atom index) of a
        residue specfied by residue id or for all atoms in a molecule (resid: all)
        '''
        try:
            assert(resid != None) #Check if resid is passed
        except AssertionError:
            sys.exit('** No residue id provided **. \nValid input to resid: residue id or all.')
            
        if str(resid).lower() == 'all':
            resid = sorted(self.residue)
        else:
            resid = [resid]
        try:
            coordinates = [self.residue[key]['atoms'][atmidx]['cord'] for key in resid for atmidx in sorted(self.residue[key]['atoms'])]
        except KeyError:
            sys.exit("*** INVALID INPUT ***. \nInput residue id out of range")
        return numpy.array(coordinates)
    
    def GetCordbyIdx(self, atmidx= None):
        '''
        Returns coordinates (as numpy array) of an atom refered by its atom index (atmidx)
        '''
        try:
            assert(atmidx != None) #Check if atom index is passed
        except AssertionError:
            sys.exit('** No atom index provided **. \nValid input to atmidx: atom index.')
        
        try:
            coordinates = next(self.residue[key]['atoms'][atmidx]['cord'] for key in sorted(self.residue) for atm_id in sorted(self.residue[key]['atoms']) if atmidx==atm_id)
        except StopIteration:
            print "*** INVALID INPUT ***. \nInput atmidx out of range"
            sys.exit("Atom index range for input "+ self.molecule_type()+ " molecule is: "+ str(self.atmidx[0])+"-"+str(self.atmidx[-1]))
        return numpy.array(coordinates)
    
    def GetDist(self, atm1, atm2):
        '''
        Returns the euclidean distance between two atoms.
        atm1 and atm2 are atom indices corresponding to atom1 and atom2
        '''
        cord_atm1 = self.GetCordbyIdx(atm1)
        cord_atm2 = self.GetCordbyIdx(atm2)
        return numpy.sqrt(((cord_atm2-cord_atm1)*(cord_atm2-cord_atm1)).sum())
    
    def GetAngle(self, atm1, atm2, atm3):
        '''
        Returns the angle (in degrees) formed by three atoms centered at atom2.
        atm1, atm2, and atm3 are atom indices corresponding to atom1, atom2, and atom3
        '''
        cord_atm1 = self.GetCordbyIdx(atm1)
        cord_atm2 = self.GetCordbyIdx(atm2)
        cord_atm3 = self.GetCordbyIdx(atm3)
        vec_21 = cord_atm1 - cord_atm2 #Vector tail at atm2 and head at atm1
        vec_23 = cord_atm3 - cord_atm2
        dot = numpy.dot(vec_21, vec_23)
        vec_21mod = numpy.sqrt((vec_21*vec_21).sum()) #length of vector 21
        vec_23mod = numpy.sqrt((vec_23*vec_23).sum())
        cos_angle = dot / vec_21mod / vec_23mod       # In radians
        return numpy.degrees(numpy.arccos(cos_angle))
    
    def GetDihed(self, atm1, atm2, atm3, atm4):
        '''
        Returns the dihedral angle (in degrees) formed by four atoms around bond connecting atm2 and atm3.
        '''
        cord_atm1 = self.GetCordbyIdx(atm1)
        cord_atm2 = self.GetCordbyIdx(atm2)
        cord_atm3 = self.GetCordbyIdx(atm3)
        cord_atm4 = self.GetCordbyIdx(atm4)
        vec_12 = cord_atm2 - cord_atm1 #Vector tail at atm1 and head at atm2
        vec_23 = cord_atm3 - cord_atm2
        vec_34 = cord_atm4 - cord_atm3
        
        uni_vec_23 = vec_23 / numpy.sqrt((vec_23*vec_23).sum()) #Unit-vector 23
        
        vec_123 = numpy.cross(vec_12, vec_23)   #Vector Normal to plane defined by atm1, atm2, and atm3
        uni_vec_123 = vec_123 / numpy.sqrt((vec_123*vec_123).sum()) #Unit-vector
        
        vec_234 = numpy.cross(vec_23, vec_34)   #Vector Normal to plane defined by atm2, atm3, and atm4
        uni_vec_234 = vec_234 / numpy.sqrt((vec_234*vec_234).sum()) #Unit-vector
        
        cross_normal = numpy.cross(uni_vec_123, uni_vec_234)
        cos_term = numpy.dot(uni_vec_123, uni_vec_234)  # In radians
        sin_term = numpy.dot(uni_vec_23, cross_normal)
        return numpy.degrees(numpy.arctan2(sin_term,cos_term))
        
        
    def GetSubset(self, criteria=''):
        '''
        Returns a subset of molecule object acccording to a given criteria.
        '''
        valid_criteria= ['ca', 'cb' ,'heavy', 'sidechain', 'backbone', 'polar', 'carbons'] #MUST be in lowercase
        protein_specific= ['ca','cb','sidechain', 'backbone'] #MUST be in lowercase
        
        #Check if input criteria is recognized
        try:
            assert criteria.lower() in valid_criteria
        except AssertionError:
            sys.exit('*** UNRECOGNIZED CRITERIA ('+str(criteria)+') SPECIFIED ***\n No Subset returned.'+
                     '\nList of recognized criteria: '+str(valid_criteria))
        
        #Check if criteria is appropriate for molecule type
        try:
            if criteria.lower() in protein_specific:
                assert self.molecule_type().lower() == 'protein'
        except AssertionError:
            sys.exit('*** PROTEIN SPECIFIC CRITERIA ('+str(criteria)+') SPECIFIED for molecule type '+str(self.molecule_type())+'***'+
                     '\n No Subset returned.')
        
        #Set select to user-specified criteria
        if criteria.lower() == 'ca':
            select = ['ca']
        elif criteria.lower() == 'cb':
            select = ['cb']
        elif criteria.lower() == 'polar':
            select = ['n', 'o']
        elif criteria.lower() in ['backbone', 'sidechain']:
            select = ['n', 'c', 'ca', 'o']
        
        OBsubset = deepcopy(self)
        #Create Subset by deleting values not required from OBsubset
        for key in sorted(OBsubset.residue):
            for atm_id in sorted(OBsubset.residue[key]['atoms']):
                if criteria.lower() not in ['sidechain', 'heavy', 'polar', 'carbons'] and OBsubset.residue[key]['atoms'][atm_id]['name'].lower() not in select:
                    del OBsubset.residue[key]['atoms'][atm_id]
                    OBsubset.atmidx.remove(atm_id)
                elif criteria.lower() == 'sidechain':
                    if OBsubset.residue[key]['atoms'][atm_id]['name'].lower() in select or OBsubset.residue[key]['atoms'][atm_id]['name'][0].lower() == 'h':
                        del OBsubset.residue[key]['atoms'][atm_id]
                        OBsubset.atmidx.remove(atm_id)
                elif criteria.lower() == 'heavy' and OBsubset.residue[key]['atoms'][atm_id]['name'][0].lower() == 'h':
                    del OBsubset.residue[key]['atoms'][atm_id]
                    OBsubset.atmidx.remove(atm_id)
                elif criteria.lower() == 'polar' and OBsubset.residue[key]['atoms'][atm_id]['name'][0].lower() not in ['n', 'o']:
                    del OBsubset.residue[key]['atoms'][atm_id]
                    OBsubset.atmidx.remove(atm_id)
                elif criteria.lower() == 'carbons' and OBsubset.residue[key]['atoms'][atm_id]['name'][0].lower() != 'c':
                    del OBsubset.residue[key]['atoms'][atm_id]
                    OBsubset.atmidx.remove(atm_id)
            
            #Check if no atoms remain then delete residue
            if not OBsubset.residue[key]['atoms']: #Empty dict evaluates to False
                del OBsubset.residue[key]
        
        #if OBsubset.molecule_type().lower() != 'ligand':
        OBsubset.resids = sorted(OBsubset.residue)
        OBsubset.nor = len(OBsubset.resids)
        #Return OBsubset only if any atom satisfying the criteria was found
        if not OBsubset.residue:
            print "NO ATOM FOUND with input criteria: %s" % criteria
        else:
            return OBsubset
        
    def GetRange(self):
        '''
        Returns min and max values of atom coordinates in the x, y, and z directions
        '''
        cord = self.GetCord(resid= 'all')
        return numpy.amin(cord, axis=0), numpy.amax(cord, axis=0)
    
    def GetMass(self):
        '''
        Returns a numpy array of mass of all atoms in the molecule
        '''
        mass_file = open(os.path.join(os.path.dirname(__file__), 'mass'), 'r')
        atom_mass={}
        for line in mass_file:
            line = line.split()
            atom_mass[line[0].upper()] = float(line[1])
        mass_file.close()
        
        try:
            return_mass = [atom_mass[self.residue[key]['atoms'][atmidx]['name']] for key in sorted(self.residue) for atmidx in sorted(self.residue[key]['atoms'])]
        except KeyError as missing:
            print '*** Missing Atom: %s' % missing+ ' in file '+ os.path.join(os.path.dirname(__file__), "mass")+' ***\n'
            sys.exit('Add missing atom to file mass and rerun!')
        return numpy.array(return_mass)
            
    def GetCenter(self, cen_type = 'geom', custom_weight=None):
        '''
        Returns the geometric center or center of mass of a molecule.
        cen_type: either geom (geometric) or mass or custom. When custom, a list of size N (no. of atoms) must be passed containing
                  weights for each atom.
        '''
        cord = self.GetCord(resid= 'all')
        
        if cen_type.lower() == 'geom':
            weight = numpy.ones(numpy.shape(cord)[0])
        elif cen_type.lower() == 'mass':
            weight = self.GetMass()
        elif cen_type.lower() == 'custom':
            try:
                assert custom_weight != None
            except AssertionError:
                sys.exit('*** Required argument custom_weight not passed ***')
            
            weight = numpy.array(custom_weight)
            try:
                assert len(weight) == numpy.shape(cord)[0]
            except AssertionError:
                sys.exit('*** Number of custom_weight passed not equals number of atoms ***')
        else:
            sys.exit('*** INVALID CENTER TYPE ***. \nArgument type must be one of [geom, center, custom].')
        cd=numpy.sum(weight)
        ab=numpy.dot(weight,cord)
        
        return numpy.dot(weight,cord)/numpy.sum(weight)
        
    
    @abstractmethod
    def molecule_type(self):
        """Return a string representing the type of molecule this is."""
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
    

    

        
    