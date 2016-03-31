#!/usr/bin/env python

from __future__ import division
import mdtraj as md
import os.path
import numpy
import sift
from Bython.Structure.configstruc import Mol_types

class Trajectory:
    '''
    Base class for MD simulation trajectories.
    Since MDTraj is used, Trajectory has access to all MDTraj attributes and functions.
    '''
    
    def __init__(self):
        self.traj_iter = None #MDTraj trajectory iterator
        self.topology = None #MDTraj topology object
        self.small_mols = () #List of small molecule ligands in trajectory
        self.name = ''
        #The varibales below are used to reload the trajectory iterator
        self.filename = None
        self.chunk = 0
        self.topfile = None #This is the name of topology file
        self.stride = 1
        self.atom_indices = None
        self.skip = 0        
                
    def Readtraj(self, filename=None, topfile=None, chunk=0, stride=1, atom_indices=None, skip=0, name=''):
        '''
        Returns MDTraj trajectory iterator
        Input:
        filename: Path to the trajectory file on disk (with file extension)
        chunk: Number of frames to load at once from disk per iteration. If 0. load all.
        top: Topology file for the trajectory (for example a .gro or .pdb file of starting structure).
        stride: Read every nth-frame, default: 1
        atom_indices: read only a subset of atom coodinates if not None.
        skip: Skipt first n frames.
        '''
        #Check if file exists
        if not os.path.isfile(filename):
            raise IOError("Cannot locate file: %s" % filename)
        if topfile != None and not os.path.isfile(topfile):
            raise IOError("Cannot locate file: %s" % topfile)
        #Call MDTraj iterload function
        self.traj_iter = md.iterload(filename=filename, chunk=chunk, top=topfile, stride=stride, atom_indices=atom_indices, skip=skip)
        #Save the read state for Reloading trajectory iterator
        self.Save_ReadState(filename=filename, chunk=chunk, topfile=topfile, stride=stride, atom_indices=atom_indices, skip=skip)
        #set the topology varibale
        for chunk in self.traj_iter:
            self.topology= chunk.topology
            break
        #Find and set small_mols if any small molecule in trajectory
        self.Sense_SmallMol()
        #Reload the instance traj_iter variable
        self.Reload()
        self.name=name
        
    def Save_ReadState(self, filename=None, chunk=0, topfile=None, stride=None, atom_indices=None, skip=0):
        '''
        Sets the instance variables to parameter values passed in Readtraj call allowing to
        reaload trajectory iterator (self.traj_iter) when it is empty.
        '''
        self.filename = filename
        self.chunk = chunk
        self.topfile = topfile
        self.stride = stride
        self.atom_indices = atom_indices
        self.skip = skip
        
    def Reload(self, newchunk=None):
        '''
        Reload the empty traj_iter instance variable to the saved read state.
        '''
        if newchunk != None:
            self.traj_iter = md.iterload(filename=self.filename, chunk=newchunk, top=self.topfile,
                                         stride=self.stride, atom_indices=self.atom_indices, skip=self.skip)
        else:
            self.traj_iter = md.iterload(filename=self.filename, chunk=self.chunk, top=self.topfile,
                                         stride=self.stride, atom_indices=self.atom_indices, skip=self.skip)
        
    def Sense_SmallMol(self):
        '''
        Sets the instance small_mols tuple to small molecule ligands in the trajectory
        '''
        self.small_mols = tuple([residue.name for residue in self.topology.residues if residue.name.lower()
                                 not in Mol_types['protein'] + Mol_types['lipid'] + Mol_types['ion'] + ('hoh',)])
        
    def SIFT(self, lig_top=None, res_index=None, bit_res=7, aro_cut=4.0, apolar_cut=4.5, hbond_cut=3.5,elec_cut=4.0,verbose=True, sc_only=False):
        '''
        For SIFT caluclation, topology must be read using .pdb format as MDTraj doesn't provide
        # bonds information for .gro file
        '''
        final_sift = [] #Final fingerprint for entire trajectory
        #Reload Trajectory in case user already worked with few/all chunks, also work with each frame at a time
        self.Reload(newchunk=1)
        loop_count = 0
        if verbose:
            print ("Structural Interaction fingerprint calculation for %s trajectory with %s ligand(s) (%s)" %
                   (self.name, str(len(self.small_mols)), ', '.join(self.small_mols)))
        for partial_traj in self.traj_iter:
            if loop_count:
                verbose=False
            #Call to appropriate bit function
            if bit_res == 7:
                partial_sift = sift.SIFT_7bit(traj_ob=partial_traj, lig_top=lig_top, res_index=res_index, aro_cut=aro_cut,
                                                   apolar_cut=apolar_cut, hbond_cut=hbond_cut,
                                                   elec_cut=elec_cut,verbose=verbose, use_sc_only=sc_only)
                final_sift.append(partial_sift)
                loop_count+=1
        print numpy.array(final_sift).shape
                
                
                
            
        
        
        
        
        
        
        
        
        
        
        
        