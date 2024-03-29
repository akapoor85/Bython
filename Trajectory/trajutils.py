#!/usr/bin/env python

from __future__ import division
import numpy
import os
from sklearn.cluster import DBSCAN, AgglomerativeClustering

def RMSD(traj, ref, idx):
    '''
    Returns RMSD between trajectory frames in mdtraj traj and reference frame in mdtraj ref object.
    RMSD is computed for atom indices specified in idx.
    '''
    return numpy.sqrt(numpy.sum(numpy.square(traj.xyz[:,idx,:] - ref.xyz[:,idx,:]),axis=(1, 2))/len(idx))

def Cluster_Traj(distance_mat=None, eps=None, min_samples=None, metric= 'precomputed', use_algo='DBSCAN',
                 n_clusters=2, linkage='average' ):
    '''
    Clustering of input Distance matrix, for ex., RMSD_Matrix
    '''
    if use_algo.lower() == 'dbscan':
        try:
            assert eps != None and min_samples != None
        except AssertionError:
            raise ValueError('DBSCAN parameters eps and min_samples cannot be None.')
        cluster_estimator = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(distance_mat)
    elif use_algo.lower() == 'hierarchical':
        cluster_estimator = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage, affinity=metric).fit(distance_mat)
    
    cluster_labels = cluster_estimator.labels_
    return cluster_labels

def Process_Clusters(traj, cluster_labels, align_index=None, rmsd_index=None, save_nFrames=100, beta=1.0, rmsd_unit_factor=10,save_format='xtc',save_path=''):
    '''
    Characterizes each cluster and saves sav_nFrames number of frames as separate traj file.
    Input: traj-> mdtraj traj object
           cluster_labels -> corresponding to each trajectory frame
           dist_mat -> Distance matrix used in clustering
           save_nFrames -> Number of trajectory frames o save for each cluster
    '''
    try:
        assert align_index != None
    except AssertionError:
        raise ValueError('align_index cannot be None')
    #If rmsd_index not defined, assign align_index to rmsd_index
    if rmsd_index == None:
        rmsd_index = align_index
    print '-'*110
    print('Total frames: %d\t Total Noise: %d\t Total Classified: %d\t # of Clusters: %d\n'
          % (len(cluster_labels), len(cluster_labels[cluster_labels == -1]), len(cluster_labels)-len(cluster_labels[cluster_labels == -1]),
             len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)))
    print "Cluster_Label\tCluster_Size\tMax_RMSD_from_Centroid\tCluster_Spread (Avg. RMSD)  Centroid_Frame_Id"
    print '-'*110
    for c_label in set(cluster_labels):
        if c_label != -1: #Ignore noise
            frame_indices = numpy.where(cluster_labels == c_label)[0]
            #All Frames corresponding to this cluster from original trajectory
            traj_cluster_frames = traj[frame_indices]
            #Compute All pairwise RMSD required for cluster medoid calculation
            
            rmsd_mat = numpy.empty((traj_cluster_frames.n_frames, traj_cluster_frames.n_frames))
            for f_index in range(traj_cluster_frames.n_frames):
                traj_cluster_frames.superpose(traj_cluster_frames[f_index], atom_indices = align_index)
                rmsd_mat[f_index]= RMSD(traj_cluster_frames, traj_cluster_frames[f_index], rmsd_index)*rmsd_unit_factor
            
            #Compute Similarity Score and return max(Similarity score) as medoid
            medoid_index = numpy.exp(-beta*rmsd_mat / rmsd_mat.std()).sum(axis=1).argmax()
            #medoid_index = 0 
            c_medoid = traj_cluster_frames[medoid_index]
            c_spread = numpy.mean(rmsd_mat[medoid_index])
            max_rmsd = numpy.max(rmsd_mat[medoid_index])
            #Save Cluster centroid pdb
            c_medoid.save(os.path.join(save_path, 'Cluster_'+str(c_label)+'_Centroid.pdb'))
            print '%d\t\t%d\t\t%f\t\t%f\t\t\t%d' % (c_label, len(frame_indices), max_rmsd,c_spread, frame_indices[medoid_index]) 
            #Superpose savenFrames cluster frames on centroid and save savenFrames cluster frames
            if save_nFrames < len(frame_indices):
                frame_indices = numpy.random.choice(frame_indices, size=save_nFrames, replace = False)
            #Superpose cluster frames on Cluster centroid; Save Cluster centroid and cluster frames
            save_cluster_frames = traj[frame_indices]
            save_cluster_frames.superpose(c_medoid, atom_indices = align_index)
            save_cluster_frames.save(os.path.join(save_path, 'Cluster_'+str(c_label)+'.'+save_format))
            
def cmdscale(distance_matrix):
    '''                                                                                       
    Classical multidimensional scaling (MDS)                                                                                       
    Input-> distance_matrix : (N, N) array                                                                          
    Output-> Y : (N, P) array                                                                          
                Configuration matrix. Each column represents a dimension. Only the                    
                p dimensions corresponding to positive eigenvalues of B are returned.                 
                Note that each dimension is only determined up to an overall sign,                    
                corresponding to a reflection.                                                        
             e : (n,) array
                 Eigenvalues of B.
    '''
    # Number of points                                                                        
    n = len(distance_matrix) 
    # Centering matrix                                                                        
    H = numpy.eye(n) - numpy.ones((n, n))/n 
    # YY^T                                                                                    
    B = -H.dot(distance_matrix**2).dot(H)/2 
    # Diagonalize                                                                             
    evals, evecs = numpy.linalg.eigh(B)
    # Sort by eigenvalue in descending order                                                  
    idx   = numpy.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]
    # Compute the coordinates using positive-eigenvalued components only                      
    w, = numpy.where(evals > 0)
    L  = numpy.diag(numpy.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L) 
    return Y, evals

def tanimoto_dissimilarity(fingerprint_matrix, relative=None):
    '''
    Borrowed from Davide's HREX analysis code, modified a bit to cover the case
    of all 0 fingerprints (commented out). Tanimoto distance between two fingerprints with all 0 bits
    is set to 0. 
    '''
    featuresum= numpy.sum(fingerprint_matrix, axis=1) # Number of on bits in each frame
    '''
    #Don't remember the logic behind this step but adding one doesn't seem right, so commenting it out
    #Check if any of the fingerprint is all zeros, and set varibales accordingly
    if not numpy.count_nonzero(featuresum) == featuresum.shape[0]: #Evaluates to True if there is a fingerprint with all 0 bits
        eps_regularize = 0 
        add_one = 1
    else:
        eps_regularize = 1e-5 
        add_one = 0
    '''
    add_one = 0
    eps_regularize= 1e-5
    ab= numpy.dot(fingerprint_matrix, fingerprint_matrix.T)
    a=numpy.repeat([featuresum],len(featuresum),axis=0)
    b=a.T
    r=1-(ab+add_one)/(a+b+add_one-ab+eps_regularize)
    if relative:
        s=r/(1-r+relative)
    else:
        s=r
    if eps_regularize:
        numpy.fill_diagonal(s, 0.0)
    return s        
    
    
    
    
    
