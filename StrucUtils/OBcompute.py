#!/usr/bin/env python
'''
Implements functions to perform calculations on Molecule_type objects
'''

from __future__ import division
import numpy
import measure
import sys

def OBDistMat(mol_object):
    '''
    Returns a N X N distance matrix where N is the number of atoms in input molecule type object
    '''
    try:
        assert mol_object.molecule_type() in ['Protein', 'Lipid', 'Ligand']
    except AttributeError:
        sys.exit("***INVALID INPUT TYPE***.\nOBDistMat expects a molecule_type (Protein, Lipid, Ligand) object.")
    
    return measure.DistMat(mol_object.GetCord('all'))    

def OBRotate(mol_object=None, axis=None, thetha=None):
    '''
    This function applies rotation matrix on input object along the axis specified in 'axis'
    and returns the new set of coordnates as NX3 array.
    'axis' can be either of x, y, or z in which case one single rotation is performed, or it
    can be a list of axis (for ex: [z,x,y,x...]) in which case multiple rotations are performed
    in the order of axis values in 'axis'.
    'thetha' is the angle of rotation in degrees, It can be a single value or a list of values only if axis
    is also a list    
    '''
    try:
        assert mol_object.molecule_type() in ['Protein', 'Lipid', 'Ligand']
    except AttributeError:
        sys.exit("***INVALID INPUT TYPE***.\nFunction expects a molecule_type (Protein, Lipid, Ligand) object.")
    
    try:
        assert(axis != None and thetha !=None) #Check if axis is passed
    except AssertionError:
        sys.exit('*** AXIS and/or THETHA CANNOT BE NONE ***. \nValid input to axis is one of: x, y, z or [list of axis in the order to perform multiple rotations].')
    
    return measure.Rotate(mol_object.GetCord('all'), axis, thetha)

"""
def OBContacts(**kwargs):
    '''
    Calculates all contacts within a molecule (if only one molecule passed)
    or between two molecules (if two molecule objects passed) and optionally writes output to a file.
    Required Input: object1= pass_a_molecule_object
    Optional Input: object2= pass_another_molecule_object
                    output= Name_of_file
    '''
    try:
        object1 = kwargs.pop('object1')
        object2 = kwargs.pop('object2', None)
        output =  kwargs.pop('output', None)
        assert len(kwargs) == 0
    except KeyError as missing:
        sys.exit("*** Required Argument Name: "+ str(missing) +" Not Passed ***" +
                 "\nVALID INPUT FORMAT: object1= molecule_object, object2= molecule_object, output= file_name")
    except AssertionError:
        sys.exit("Unrecognized argument name passed: %s" % ",".join(kwargs.keys()))
"""     
    
        
        
            
