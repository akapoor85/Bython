#!/usr/bin/env python

from __future__ import division
import numpy
import sys


def DistMat(cord):
    '''
    Equivalent to Matlab pdsit function. It takes input a NX3 numpy array (cord)
    and returns a NXN matrix where the ijth cell is the distance between the ith
    and jth point in NX3 matrix
    '''
    d0 = numpy.subtract.outer(cord[:,0], cord[:,0])**2
    d1 = numpy.subtract.outer(cord[:,1], cord[:,1])**2
    d2 = numpy.subtract.outer(cord[:,2], cord[:,2])**2
    dist = numpy.sqrt(d0+d1+d2)
    return dist

def Rotate(cord, axis=None, thetha=None):
    '''
    This function applies rotation matrix along the axis specified in 'axis'
    on the input NX3 numpy array (cord) and returns the new set of coordnates as NX3 array.
    'axis' can be either of x, y, or z in which case one single rotation is performed, or it
    can be a list of axis (for ex: [z,x,y,x...]) in which case multiple rotations are performed
    in the order of axis values in 'axis'.
    'thetha' is the angle of rotation in degrees, It can be a single value or a list of values only if axis
    is also a list    
    '''
    try:
        assert(axis != None and thetha !=None) #Check if axis is passed
    except AssertionError:
        sys.exit('*** AXIS and/or THETHA CANNOT BE NONE ***. \nValid input to axis is one of: x, y, z or [list of axis in the order to perform multiple rotations].')
    
    if not type(axis) is list:
        axis = [axis]
    if not type(thetha) is list:
        thetha = [thetha]    
    
    try:
        assert(len(axis) == len(thetha))
    except AssertionError:
        sys.exit('*** Number of input axis does not match number of thetha values ***')
        
    cos_term = numpy.cos(numpy.radians(numpy.array(thetha)))
    sine_term = numpy.sin(numpy.radians(numpy.array(thetha)))
    newcord = numpy.transpose(cord)
    
    for rot_axis, cos_thetha, sin_thetha in zip(axis,cos_term,sine_term):
        if rot_axis.lower() == 'x':
            rot_mat = numpy.array([[1, 0, 0], [0, cos_thetha, -sin_thetha], [0, sin_thetha, cos_thetha]])
        elif rot_axis.lower() == 'y':
            rot_mat = numpy.array([[cos_thetha, 0, sin_thetha], [0, 1, 0], [-sin_thetha, 0, cos_thetha]])
        elif rot_axis.lower() == 'z':
            rot_mat = numpy.array([[cos_thetha, -sin_thetha, 0], [sin_thetha, cos_thetha, 0], [0, 0, 1]])
        else:
            sys.exit('*** INVALID AXIS ARGUMENT %s ***' % rot_axis+'. \nValid input to axis is one of: x, y, z or [list of x,y,z ].')
        
        newcord = numpy.dot(rot_mat, newcord)
    return numpy.transpose(newcord)
            
    
