#!/usr/bin/env python

from __future__ import division
import mdtraj as md
import os.path, sys
import numpy
import time
import sift
import sift_pro_pro
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
        
    def Reload(self, newchunk=None, newskip=0, newstride=1):
        '''
        Reload the empty traj_iter instance variable to the saved read state.
        '''
        if newchunk != None:
            self.chunk = newchunk
        if newskip != 0:
            self.skip = newskip
        if newstride !=1:
            self.stride = newstride
        
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
                    
    def SIFT(self, lig_top=None, res_index=None, bit_res=7, aro_cut=4.0, apolar_cut=4.5,
             hbond_cut=3.5,elec_cut=4.0,verbose=True, sc_only=False, ligand_neighb_dist=None, res_index_include_all=[], add_bonds=[]):
        '''
        Ligand-protein Structural Interaction Fingerprint calculation.
        Fingerprint can be calculated for any number of ligands together. Currently supports 7-bit and 9-bit resolution per protein residue.
        SIFT can be calculated with all residues (or specific residues) and all-atoms (or just the protein sidechain) for the
        entire trajectory.
        For SIFT caluclation, topology must be read using .pdb format as MDTraj doesn't provide bonds information for .gro file
        res_index: Is a 0-based index i.e. res_index = Residue id-1
        ligand_neighb_dist: Residues within 'ligand_neighb_dist' Ang of ligand polar atoms will be used for water-mediated hydrogen bond determination
                       (default: (3*hbond_cut+0.05) Ang)
	res_index_include_all: For these residue indexes include both sidechain and backbone atoms even if sc_only is True. 
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
        if ligand_neighb_dist == None:
            ligand_neighb_dist = 3*hbond_cut + 0.05
        if verbose:
            print ("Structural Interaction fingerprint calculation for %s trajectory with %s ligand(s) (%s)" %
                   (self.name, str(len(self.small_mols)), ', '.join(self.small_mols)))
        s_time = time.clock()
        
	for partial_traj in self.traj_iter:
	    #Add bonds to topology if needed
	    if len(add_bonds):
		self.Add_Bonds_to_Topology(add_bonds, mdtraj_object=partial_traj)

            if loop_count:
                verbose=False
            #Call to appropriate bit function
            if bit_res == 7:
                partial_sift, protein_resindex = sift.SIFT_7bit(traj_ob=partial_traj, lig_top=lig_top, res_index=res_index, aro_cut=aro_cut,
                                                   apolar_cut=apolar_cut, hbond_cut=hbond_cut,
                                                   elec_cut=elec_cut,verbose=verbose, use_sc_only=sc_only, res_index_include_all=res_index_include_all)
            elif bit_res == 9:
                partial_sift, protein_resindex = sift.SIFT_9bit(traj_ob=partial_traj, lig_top=lig_top, res_index=res_index, aro_cut=aro_cut,
                                                   apolar_cut=apolar_cut, hbond_cut=hbond_cut,
                                                   elec_cut=elec_cut,verbose=verbose, use_sc_only=sc_only, ligand_neighb_dist=ligand_neighb_dist, res_index_include_all=res_index_include_all)
            
            final_sift.append(partial_sift)
            loop_count+=1
                      
        e_time = time.clock()
        print "SIFT Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.2f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/numpy.array(final_sift).shape[0])
        #Prepare protein residue labels
        residue_name_num = [partial_traj.topology.residue(residue_index) for residue_index in protein_resindex]
        residue_labels = []
        
        for each_res in residue_name_num:
            if sc_only and str(each_res)[0:3].upper() == 'GLY' and not each_res.index in res_index_include_all:
                continue
            elif bit_res==7:
                for sift_type in ['Apolar', 'Aro_F2F', 'Aro_E2F', 'Hbond_ProD', 'Hbond_ProA', 'Elec_ProP', 'Elec_ProN']:
                    residue_labels.append(str(each_res)+'_'+sift_type)
            elif bit_res==9:
                for sift_type in ['Apolar', 'Aro_F2F', 'Aro_E2F', 'Hbond_ProD', 'Hbond_ProA', 'Elec_ProP', 'Elec_ProN', 'Hbond_1Wat', 'Hbond_2Wat']:
                    residue_labels.append(str(each_res)+'_'+sift_type)
        #Calculate Tanimoto dissimilarity matrix
        self.tanimoto_dissimilarity_matrix = trajutils.tanimoto_dissimilarity(numpy.array(final_sift))
        
        return numpy.array(final_sift), residue_labels, numpy.mean(numpy.array(final_sift), axis=0)
    
    def SIFT_PROPRO(self, res_index=None, res_index2=None, bit_res=7, aro_cut=4.0, apolar_cut=4.5,
             hbond_cut=3.5,elec_cut=4.0,verbose=True, sc_only=False):
        '''
        protein-protein Structural Interaction Fingerprint calculation.
        Fingerprint can be calculated for any number of ligands together. Currently supports 7-bit and 9-bit resolution per protein residue.
        SIFT can be calculated with all residues (or specific residues) and all-atoms (or just the protein sidechain) for the
        entire trajectory.
        For SIFT caluclation, topology must be read using .pdb format as MDTraj doesn't provide bonds information for .gro file
        res_index: Is a 0-based index i.e. res_index = Residue id-1
        ligand_neighb_dist: Residues within 'ligand_neighb_dist' Ang of ligand polar atoms will be used for water-mediated hydrogen bond determination
                       (default: (3*hbond_cut+0.05) Ang)
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
    
        s_time = time.clock()
        for partial_traj in self.traj_iter:
            if loop_count:
                verbose=False
            #Call to appropriate bit function
            if bit_res == 7:
                partial_sift, protein_resindex = sift_pro_pro.SIFT_7bit(traj_ob=partial_traj, res_index2=res_index2, res_index=res_index, aro_cut=aro_cut,
                                                   apolar_cut=apolar_cut, hbond_cut=hbond_cut,
                                                   elec_cut=elec_cut,verbose=verbose, use_sc_only=sc_only)
            elif bit_res == 9:
                partial_sift, protein_resindex = sift_pro_pro.SIFT_9bit(traj_ob=partial_traj, lig_top=lig_top, res_index=res_index, aro_cut=aro_cut,
                                                   apolar_cut=apolar_cut, hbond_cut=hbond_cut,
                                                   elec_cut=elec_cut,verbose=verbose, use_sc_only=sc_only, ligand_neighb_dist=ligand_neighb_dist)
            
            final_sift.append(partial_sift)
            loop_count+=1
                      
        e_time = time.clock()
        print "SIFT Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.2f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/numpy.array(final_sift).shape[0])
        #Prepare protein residue labels
        residue_name_num = [partial_traj.topology.residue(residue_index) for residue_index in protein_resindex]
        residue_labels = []
        
        for each_res in residue_name_num:
            if sc_only and str(each_res)[0:3].upper() == 'GLY':
                continue
            elif bit_res==7:
                for sift_type in ['Apolar', 'Aro_F2F', 'Aro_E2F', 'Hbond_ProD', 'Hbond_ProA', 'Elec_ProP', 'Elec_ProN']:
                    residue_labels.append(str(each_res)+'_'+sift_type)
            elif bit_res==9:
                for sift_type in ['Apolar', 'Aro_F2F', 'Aro_E2F', 'Hbond_ProD', 'Hbond_ProA', 'Elec_ProP', 'Elec_ProN', 'Hbond_1Wat', 'Hbond_2Wat']:
                    residue_labels.append(str(each_res)+'_'+sift_type)
        #Calculate Tanimoto dissimilarity matrix
        self.tanimoto_dissimilarity_matrix = trajutils.tanimoto_dissimilarity(numpy.array(final_sift))
        
        return numpy.array(final_sift), residue_labels, numpy.mean(numpy.array(final_sift), axis=0)
    
    def RMSD_Matrix(self, align_index=None, rmsd_index=None, ref_frame=None, ref_md_ob=None,rmsd_unit_factor=10,pre_aligned=False):
        '''
        Sets the class instance variable rmsd_matrix (NxN array) containing RMSD (in Angstrom) between each of the N trajectory frames.
        If ref_frame is not none, then a 1xN rmsd array is set with rmsd of all frames w.r.t ref_frame of trajectory in self or
        the ref_frame of a mdtraj reference object (if provided as input)  
        !!Note: Superpose function will modify the original coordinates in traj.
        To prevent this, pass a deepcopy image of the object.!! 
        Input: align_index-> list of atom indices used in structure alignment (0-based).
               rmsd_index-> list of atom indices used for calculating rmsd (0-based).
                            If not specified, align_index will be used to compute rmsd
               ref_frame-> Specific reference frame index to use. Frames number from 0 to N-1.
               ref_md_ob-> Compute rmsd w.r.t ref_frame of trajectory reffered by ref_md_ob
               rmsd_unit_factor: Default unit in mdtraj is nanometers, here RMSD is returned in Angstrom.
               pre_aligned -> skips superpose step if trajectory was aligned already. Default is to align the trajectory first.
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
                if ref_md_ob != None:
                    if not pre_aligned:
                        full_traj.superpose(ref_md_ob[ref_frame],  atom_indices = align_index)
                    rmsd_mat = trajutils.RMSD(full_traj, ref_md_ob[ref_frame], rmsd_index)*rmsd_unit_factor
                else:
                    if not pre_aligned:
                        full_traj.superpose(full_traj[ref_frame],  atom_indices = align_index)
                    rmsd_mat = trajutils.RMSD(full_traj, full_traj[ref_frame], rmsd_index)*rmsd_unit_factor
            else:
                #Otherwise compute and return pairwise RMSDs
                rmsd_mat = numpy.empty((full_traj.n_frames, full_traj.n_frames))
                for frame_index in range(full_traj.n_frames):
                    if not pre_aligned:
                        full_traj.superpose(full_traj[frame_index],  atom_indices = align_index)
                    rmsd_mat[frame_index]= trajutils.RMSD(full_traj, full_traj[frame_index], rmsd_index)*rmsd_unit_factor
        e_time = time.clock()
        print "RMSD Matrix Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.2f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/full_traj.n_frames)
        self.rmsd_matrix = rmsd_mat
    
    def Cluster_Trajectory(self, eps=None, min_samples=None, use_algo='DBSCAN', metric='precomputed',
                           use_dist='rmsd', n_clusters=2, linkage='average', custom_dist=None):
        '''
        Clusters trajectory using input algorithm (currently only DBSCAN) and input distance matrix format (rmsd or tanimoto)
        '''
        valid_algo = ['DBSCAN', 'HIERARCHICAL']
        valid_dist = ['rmsd', 'tanimoto', 'custom']
        try:
            assert use_algo.upper() in valid_algo and use_dist.lower() in valid_dist
        except AssertionError:
            raise ValueError('Valid values for parameters\nuse_algo: %s\nuse_dist: %s' % (', '.join(valid_algo), ', '.join(valid_dist)))
        
        if use_dist.lower() == 'rmsd':
            distance_mat = self.rmsd_matrix
        elif use_dist.lower() == 'tanimoto':
            distance_mat = self.tanimoto_dissimilarity_matrix
        elif use_dist.lower() == 'custom':
            distance_mat = custom_dist
        
        self.traj_cluster_labels = trajutils.Cluster_Traj(distance_mat=distance_mat, eps=eps,
                                                          min_samples=min_samples, metric= metric,
                                                          use_algo=use_algo, n_clusters=n_clusters, linkage=linkage)
        
                
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
            
    def Get_Distance(self, atom_pairs, periodic=True, newchunk=100, unitfactor=10):
        '''
        Uses MDTraj compute_distances to compute the distances between pairs of atoms in each frame.
        Input: atom_pairs-> An array of shape (num_pairs,2) with each row giving indices of two atoms to compute distnaces for.
               periodic-> If periodic is True and the trajectory contains unitcell information, compute distances under
               the minimum image convention. Default distance units for MDtraj is in nm. Results are returned in Angstrom.  
        Returns: an array of shape (num_frames, num_pairs) with the distances, in each frame, between each pair of atoms.
        '''
        dist_pairs = numpy.empty((0, numpy.shape(atom_pairs)[0])) #Return array
        #Reload Trajectory in case user already worked with few/all chunks, working with newchunk frames at a time
        self.Reload(newchunk=newchunk)
        
        s_time = time.clock()
        for partial_traj in self.traj_iter:
            partial_dist = md.compute_distances(partial_traj, atom_pairs, periodic=periodic)
            dist_pairs = numpy.vstack((dist_pairs,partial_dist))
        e_time = time.clock()
        print "Distance Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.4f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/numpy.array(dist_pairs).shape[0])
        
        return numpy.array(dist_pairs)*unitfactor
            
    def RMSF(self, align_index=None, rmsf_index=None, ref_frame=None, rmsd_unit_factor=10):
        '''
        Returns an array of shape (rmsf_index,) containing root mean square fluctuation (in Angstrom) of each atom
        in rmsf_index. By default, average structure is used as reference. Trajectory frame index specified by ref_frame
        is used instead when ref_frame is not None. An array containing residue id is also returned. 
        !!Note: Superpose function will modify the original coordinates in traj.
        To prevent this, pass a deepcopy image of the object.!! 
        Input: align_index-> list of atom indices used in structure alignment (0-based).
               rmsf_index-> list of atom indices used for calculating rmsf (0-based).
                            If not specified, align_index will be used to compute rmsf
               ref_frame-> Specific reference frame index to use. Frames number from 0 to N-1.
               rmsd_unit_factor: Default unit in mdtraj is nanometers, here RMSD is returned in Angstrom.
        '''
        try:
            assert align_index != None
        except AssertionError:
            raise ValueError('align_index cannot be None')
        #If rmsf_index not defined, assign align_index to rmsf_index
        if rmsf_index == None:
            rmsf_index = align_index
        #Reload Trajectory in case user already worked with few/all chunks, also work with entire trajectory at a time
        self.Reload(newchunk=0)
        
        s_time = time.clock()
        
        for full_traj in self.traj_iter: #This loop runs only once
            #Get the residue id of atoms in rmsf_index
            print full_traj.n_frames
            residue_ids = [full_traj.topology.atom(atom_idx).residue.resSeq for atom_idx in rmsf_index]
            #Align trajecotry on frame 0; rmsf calculation is independent of which frame traj is aligned to.
            full_traj.superpose(full_traj[0],  atom_indices=align_index)
            #If reference frame specified, then return RMSF with the reference frame.
            if ref_frame != None:
                ref_pos = full_traj.xyz[ref_frame,rmsf_index,:]
            else:
                #Compute average structure as reference frame
                ref_pos = full_traj.xyz[:,rmsf_index,:].mean(0)
            #Compute RMSF
            diff = full_traj.xyz[:,rmsf_index,:] - ref_pos
            root_mean_sq_diff = (((diff**2).sum(2)).mean(0))**0.5         
            
        e_time = time.clock()
        print "RMSF Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.2f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/full_traj.n_frames)
        return root_mean_sq_diff*rmsd_unit_factor, residue_ids
    
    def Get_Dihedral(self, atom_indices, periodic=True, newchunk=100, unit='degrees'):
        '''
        Uses MDTraj compute_dihedrals to compute the dihedrals formed by input atom indices in each frame.
        Input: atom_indices-> An array of shape (num_dihedrals,4) with each row giving indices of four atoms to compute dihedrals for.
                              The angle is between the planes spanned by the first 3 atoms and the last 3 atoms, a torsion around the
                              bond between the middle two atoms. 
               periodic-> If periodic is True and the trajectory contains unitcell information, minimum image convention will be used
                          for dihedrls crossing periodic boundary.
               unit -> MDTraj default unit is radian. Here by default, angles are returned in degrees.
        Returns: an array of shape (num_frames, num_dihedrals).
        '''
        dihed_angles = numpy.empty((0, numpy.shape(atom_indices)[0])) #Return array
        #Reload Trajectory in case user already worked with few/all chunks, working with newchunk frames at a time
        self.Reload(newchunk=newchunk)
        
        s_time = time.clock()
        for partial_traj in self.traj_iter:
            partial_dihed = md.compute_dihedrals(partial_traj, atom_indices, periodic=periodic)
            dihed_angles = numpy.vstack((dihed_angles,partial_dihed))
        e_time = time.clock()
        print "Dihedral Calculation Finished! \nTotal time spent: \n\t%0.2f seconds process time\
        \n\tOn average, %0.4f seconds per trajectory frame" % (e_time-s_time, (e_time-s_time)/numpy.array(dihed_angles).shape[0])
        
        if unit.lower() == 'degrees':
            return numpy.degrees(dihed_angles)
        else:
            return dihed_angles

    def Get_Atom_Object(self, res_name = None, atom_name = None):
	target_residue_found = False
	for target_residue in self.topology.residues:
		if res_name == (target_residue.name+str(target_residue.index)).lower():
			target_residue_found = True
			#save atom object
			target_atom_found = False
			for target_atom in target_residue.atoms:
				if atom_name == target_atom.name.lower():
					target_atom_found = True
					atom_ob = target_atom
					break
			if not target_atom_found:
				raise ValueError('Atom %s was not found in residue %s' % (atom_name, res_name))
	if not target_residue_found:
		raise ValueError('Residue %s was not found in topology' % (res_name))
	return atom_ob
    
    def Add_Bonds_to_Topology(self, bond_array=None, mdtraj_object=None):
        '''
        Add new bonds to topology file.
        Input: bond_array-> An array of shape (num_of_bonds,2) with each row containing information on two atoms involved in bond.
                            Each column contains atom specification in the format ResnameIndex-Atomname, for ex: TYR1-N.
			    Note the residue index is the MDtraj resindex property which starts from 0 and is unique for each residue unlike PDB residue number.
        '''
	#Add each bond one at a time
	for each_bond in bond_array:
		first_res, first_res_atom = each_bond[0].split('-')[0].lower(), each_bond[0].split('-')[1].lower()
		second_res, second_res_atom = each_bond[1].split('-')[0].lower(), each_bond[1].split('-')[1].lower()
		#Get atom objects of first and second residue
		atom1_ob = self.Get_Atom_Object(first_res, first_res_atom)
		atom2_ob = self.Get_Atom_Object(second_res, second_res_atom)
		#Add bond to topology
		self.topology.add_bond(atom1_ob, atom2_ob)
		if mdtraj_object !=None:
			mdtraj_object.topology.add_bond(atom1_ob, atom2_ob)

        
            
        
        
        
    
    
    
               
                    
                
            
        
        
        
        
        
        
        
        
        
        
        
        
