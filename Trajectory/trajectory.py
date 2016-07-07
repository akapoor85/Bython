#!/usr/bin/env python

from __future__ import division
import mdtraj as md
import os.path, sys
import numpy
import time
import sift
import trajutils
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
        #Function specific variables
        self.rmsd_matrix = []
        self.tanimoto_dissimilarity_matrix = []
        self.traj_cluster_labels = []
                
    def Readtraj(self, filename=None, topfile=None, chunk=0, stride=1, atom_indices=None, skip=0, name=''):
        '''
        Returns MDTraj trajectory iterator
        Input:
        filename: Path to the trajectory file on disk (with file extension)
        chunk: Number of frames to load at once from disk per iteration. If 0, load all.
        top: Topology file for the trajectory (for example a .gro or .pdb file of starting structure).
        For SIFT caluclation, topology must be read using .pdb format as MDTraj doesn't provide bonds information for .gro file
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
    
    def Get_Index(self, atom_criteria='ca', res_index_list = []):
        '''
        Returns indices of atoms as per the atom_criteria for residues in res_index_list.
        Residue index is 0-based.
        Note: MDTraj is_sidechain does not work for hydrogen atoms, so I kept only sidechain-heavy in atom_criteria
        resn_3LETTERCODE: only returns heavy atom index
        '''
        valid_criteria= ['ca', 'cb' ,'heavy', 'sidechain-heavy', 'backbone', 'polar', 'carbons', 'resn_3LETTERCODE'] #MUST be in lowercase
        #protein_specific= ['ca','cb','sidechain', 'backbone'] #MUST be in lowercase
        #Check if input criteria is recognized
        try:
            assert atom_criteria.lower() in valid_criteria
        except AssertionError:
            if not atom_criteria.lower()[0:5] == 'resn_':
                raise ValueError('UNRECOGNIZED CRITERIA ('+str(atom_criteria)+') SPECIFIED\n No Subset returned.'+
                                 '\nList of recognized criteria: '+str(valid_criteria))
        #Set select to user-specified criteria
        if atom_criteria.lower() == 'ca':
            select = ['ca']
        elif atom_criteria.lower() == 'cb':
            select = ['cb']
        elif atom_criteria.lower() == 'polar':
            select = ['n', 'o']
        elif atom_criteria.lower() in ['backbone', 'sidechain-heavy']:
            select = ['n', 'c', 'ca', 'o']
        elif atom_criteria.lower()[0:4] == 'resn':
            select = [atom_criteria.split('_')[1].lower()]
            
        #If res_index_list is empty and specific residue index not required than use all protein residues
        if not res_index_list:
            if atom_criteria.lower()[0:4] != 'resn':
                res_index_list = [atom_object.residue.index for atom_object in self.topology.atoms if atom_object.name.lower() == 'ca']
            else:
                res_index_list = range(self.topology.n_residues)
        #Go through all-atoms in the topology
        selected_atom_indices = []
        for atom_object in self.topology.atoms:
            if atom_object.residue.index in res_index_list:
                if atom_criteria.lower() not in ['sidechain-heavy','heavy', 'polar', 'carbons'] and atom_object.name.lower() in select:
                    selected_atom_indices.append(atom_object.index)
                elif atom_criteria.lower() == 'sidechain-heavy' and (atom_object.is_sidechain == True and atom_object.element.symbol != 'H'):
                    selected_atom_indices.append(atom_object.index)
                elif atom_criteria.lower()[0:4] == 'resn':
                    if atom_object.residue.name.lower() in select and atom_object.element.symbol != 'H':
                        selected_atom_indices.append(atom_object.index)
        return selected_atom_indices
                
                    
    def SIFT(self, lig_top=None, res_index=None, bit_res=7, aro_cut=4.0, apolar_cut=4.5, hbond_cut=3.5,elec_cut=4.0,verbose=True, sc_only=False):
        '''
        Ligand-protein Structural Interaction Fingerprint calculation.
        Fingerprint can be calculated for any number of ligands together. Currently supports 7-bit resolution per protein residue.
        SIFT can be calculated with all residues (or specific residues) and all-atoms (or just the protein sidechain) for the
        entire trajectory.
        For SIFT caluclation, topology must be read using .pdb format as MDTraj doesn't provide bonds information for .gro file
        res_index: Is a 0-based index i.e. res_index = Residue id-1        
        Returns->
                final_sift: NxM array of fingerprints whre N=Number of trajecotry frames and M=bit_resolution*Number_of_residues*Number_of_ligands
                residue_labels: 1xM list containing labels corresponding to each fingeprint
                average_fpt: 1xM array of average fingerprint over entire trajectory.
                Additionally calculates tanimoto dissimilarity matrix and sets the class instance variable tanimoto_dissimilarity_matrix
        '''
        final_sift = [] #Final fingerprint for entire trajectory
        #Reload Trajectory in case user already worked with few/all chunks, also work with each frame at a time
        self.Reload(newchunk=1)
        loop_count = 0
        if verbose:
            print ("Structural Interaction fingerprint calculation for %s trajectory with %s ligand(s) (%s)" %
                   (self.name, str(len(self.small_mols)), ', '.join(self.small_mols)))
        s_time = time.clock()
        for partial_traj in self.traj_iter:
            if loop_count:
                verbose=False
            #Call to appropriate bit function
            if bit_res == 7:
                partial_sift, protein_resindex = sift.SIFT_7bit(traj_ob=partial_traj, lig_top=lig_top, res_index=res_index, aro_cut=aro_cut,
                                                   apolar_cut=apolar_cut, hbond_cut=hbond_cut,
                                                   elec_cut=elec_cut,verbose=verbose, use_sc_only=sc_only)
                final_sift.append(partial_sift)
                loop_count+=1
        #print numpy.array(final_sift).shape
        e_time = time.clock()
        print "SIFT Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.2f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/numpy.array(final_sift).shape[0])
        #Prepare protein residue labels
        residue_name_num = [partial_traj.topology.residue(residue_index) for residue_index in protein_resindex]
        residue_labels = []
        for each_res in residue_name_num:
            if sc_only and str(each_res)[0:3].upper() == 'GLY':
                continue
            for sift_type in ['Apolar', 'Aro_F2F', 'Aro_E2F', 'Hbond_ProD', 'Hbond_ProA', 'Elec_ProP', 'Elec_ProN']:
                residue_labels.append(str(each_res)+'_'+sift_type)
        #Calculate Tanimoto dissimilarity matrix
        self.tanimoto_dissimilarity_matrix = trajutils.tanimoto_dissimilarity(numpy.array(final_sift))
        return numpy.array(final_sift), residue_labels, numpy.mean(numpy.array(final_sift), axis=0)
    
    def RMSD_Matrix(self, align_index=None, rmsd_index=None, ref_frame=None, rmsd_unit_factor=10):
        '''
        Sets the class instance variable rmsd_matrix (NxN array) containing RMSD (in Angstrom) between each of the N trajectory frames.
        If ref_frame is not none, then a 1xN rmsd array is set with rmsd of all frames w.r.t ref_frame  
        !!Note: Superpose function will modify the original coordinates in traj.
        To prevent this, pass a deepcopy image of the object.!! 
        Input: align_index-> list of atom indices used in structure alignment (0-based).
               rmsd_index-> list of atom indices used for calculating rmsd (0-based).
                            If not specified, align_index will be used to compute rmsd
               ref_frame-> Specific reference frame index to use. Frames number from 0 to N-1.
               rmsd_unit_factor: Default unit in mdtraj is nanometers, here RMSD is returned in Angstrom.
        '''
        try:
            assert align_index != None
        except AssertionError:
            raise ValueError('align_index cannot be None')
        #If rmsd_index not defined, assign align_index to rmsd_index
        if rmsd_index == None:
            rmsd_index = align_index
        #Reload Trajectory in case user already worked with few/all chunks, also work with entire trajectory at a time
        self.Reload(newchunk=0)
        s_time = time.clock()
        #Compute pairwise RMSD
        for full_traj in self.traj_iter: #This loop runs only once
            #If reference frame specified, then return RMSD with the reference frame.
            if ref_frame != None:
                full_traj.superpose(full_traj[ref_frame],  atom_indices = align_index)
                rmsd_mat = trajutils.RMSD(full_traj, full_traj[ref_frame], rmsd_index)*rmsd_unit_factor
            else:
                #Otherwise compute and return pairwise RMSDs
                rmsd_mat = numpy.empty((full_traj.n_frames, full_traj.n_frames))
                for frame_index in range(full_traj.n_frames):
                    full_traj.superpose(full_traj[frame_index],  atom_indices = align_index)
                    rmsd_mat[frame_index]= trajutils.RMSD(full_traj, full_traj[frame_index], rmsd_index)*rmsd_unit_factor
        e_time = time.clock()
        print "RMSD Matrix Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.2f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/numpy.array(final_sift).shape[0])
        self.rmsd_matrix = rmsd_mat
    
    def Cluster_Trajectory(self, eps=None, min_samples=None, use_algo='DBSCAN', metric='precomputed', use_dist='rmsd'):
        '''
        Clusters trajectory using input algorithm (currently only DBSCAN) and input distance matrix format (rmsd or tanimoto)
        '''
        valid_algo = ['DBSCAN']
        valid_dist = ['rmsd', 'tanimoto']
        try:
            assert use_algo.upper() in valid_algo and use_dist.lower() in valid_dist
        except AssertionError:
            raise ValueError('Valid values for parameters\nuse_algo: %s\nuse_dist: %s' % (', '.join(valid_algo), ', '.join(valid_dist)))
        if use_algo.upper() == 'DBSCAN':
            if use_dist.lower() == 'rmsd':
                self.traj_cluster_labels = trajutils.Cluster_DBSCAN(distance_mat=self.rmsd_matrix, eps=eps,
                                                                    min_samples=min_samples, metric= metric)
            elif use_dist.lower() == 'tanimoto':
                self.traj_cluster_labels = trajutils.Cluster_DBSCAN(distance_mat=self.tanimoto_dissimilarity_matrix, eps=eps,
                                                                    min_samples=min_samples, metric= metric)
                
    def Characterize_Clusters(self, align_index=None,rmsd_index=None,save_nFrames=100, beta=1.0,rmsd_unit_factor=10,save_format='xtc',save_path=None):
        '''
        Processes each cluster, prints out cluster information, computes centroid conformation for each cluster and 
        writes it out as pdb file. Centroid is defined as the frame with maximum similarity score
        with other frames in the cluster (RMSD is scored). The similarity score is calculated as:
        Sij = exp(-beta*Dij/Dscale), where Dij = rmsd between frame i and j, and Dscale= standard deviation of D.
        Additionally, traj file for each cluster centroid is written with sav_nFrames number of frames selected randomly.
        Input: align_index-> list of atom indices used in structure alignment (0-based).
               rmsd_index-> list of atom indices used for calculating rmsd (0-based).
                            If not specified, align_index will be used to compute rmsd
               rmsd_unit_factor: Default unit in mdtraj is nanometers, here RMSD is returned in Angstrom (default: 10).
               save_nFrames: Number of frames to save for each cluster (default: 100).
               beta: Beta factor to use in computing similarity score (default: 1).
               save_format: Saves the trajectory in specified format (default: xtc).
               save_path: Path where files will be saved (default: directory from which the script was started).
        '''
        #Reload Trajectory, and work with entire trajectory at once
        self.Reload(newchunk=0)
        if save_path == None:
            save_path = os.path.abspath(os.path.dirname(sys.argv[0]))
        for full_traj in self.traj_iter: #This loop runs only once
            trajutils.Process_Clusters(full_traj, self.traj_cluster_labels, align_index=align_index,
                                       rmsd_index=rmsd_index, save_nFrames=save_nFrames, beta=beta,
                                       rmsd_unit_factor=rmsd_unit_factor, save_format=save_format, save_path=save_path)
        
            
        
        
        
    
    
    
               
                    
                
            
        
        
        
        
        
        
        
        
        
        
        
        