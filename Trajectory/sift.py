#!/usr/bin/env python

from __future__ import division
import mdtraj as md
import numpy
from Bython.Cheminfo import Bymol
from Bython.Structure.configstruc import NON_POLAR_ATOMS, PROTEIN_AROMATIC_RES, POLAR_ATOMS, PROTEIN_POSITIVE, PROTEIN_NEGATIVE
from itertools import product, chain

def lig_info(lig_ob_list=None):
    lig_names = [lig_t.name.upper() for lig_t in lig_ob_list]
    no_lig_aro = [str(lig_t.n_aromatic) for lig_t in lig_ob_list]
    no_heavy_lig = [lig_t.n_heavy for lig_t in lig_ob_list ] #number of heavy atoms in ligand
    lig_ring_indices = [lig_t.rings for lig_t in lig_ob_list]
    return lig_names, no_lig_aro, no_heavy_lig, lig_ring_indices

def SIFT_7bit(traj_ob=None, lig_top=None, res_index=None, aro_cut=4.0, apolar_cut=4.5, hbond_cut=3.0,elec_cut=4.0,verbose=True, use_sc_only=False):
    #Make a list of all protein residue indexes if res_index not defined
    if res_index == None:
        res_index = [traj_ob.topology.atom(atom_index).residue.index for atom_index in traj_ob.topology.select('protein and name CA')]
    #Go through Lig topologies
    if isinstance(lig_top, list): #If a list of lig_top object passed
        lig_names, no_lig_aro, no_heavy_lig, lig_ring_indices = lig_info(lig_top)
    elif isinstance(lig_top, Bymol): #Only one lig_top object passed
        lig_top = [lig_top] #Will also be used later in determining interactions
        lig_names, no_lig_aro, no_heavy_lig, lig_ring_indices = lig_info(lig_top)
    else:
        raise TypeError("Invalid argument type for lig_top (%s) passed. Expected a Bymol object." % (type(lig_top)))
        
    #Get the ligand heavy atom indices
    lig_heavy_indices = []
    for l_name in lig_names:
        temp=[]
        for a_idx in traj_ob.topology.select('resname ' + l_name):
            if traj_ob.topology.atom(a_idx).element.symbol != 'H':
                temp.append(a_idx)
        lig_heavy_indices.append(temp)
               
    if verbose:
        print ("Protein residues used: %s") % res_index
        print ("Ligand(s) used: %s") % ', '.join(lig_names)
        print ("Number of heavy atoms in Ligand(s) used: %s") % no_heavy_lig
        print ("Number of aromatic atoms in Ligand(s) used: %s") % no_lig_aro
        print ("Ring atom indices in Ligand(s) used: %s") % lig_ring_indices
        print ("Fingerprint resolution: 7-bits per residue (Apolar, Aromatic Face to Face, Aromatic Edge to Face, "+
               "Hbond Protein as Donor, Hbond Protein as Acceptor, Electrostatic Protein +ve, Electrostatic Protein -ve)")
    #Initialize the array, No_residue X 7 X No_of_ligands
    #sift7bit= numpy.zeros((len(res_index)*7*len(lig_names)), dtype=numpy.int) #traj_ob.xyz.shape[0],
    sift7bit = {} #Key: Ligand index, value: list of bits = 7*no_residues
    #Go through each residue and compute its interaction with each ligands
    for q_res in res_index:
        #Get the residue heavy atom indices
        if not use_sc_only:
            res_heavy_idx = [a_idx for a_idx in traj_ob.topology.select("resid " + str(q_res))
                             if traj_ob.topology.atom(a_idx).element.symbol != 'H']
        else:
            #Do only if residue is not Gly
            if traj_ob.topology.residue(q_res).name.lower() == 'gly':
                continue
            res_heavy_idx = [a_idx for a_idx in traj_ob.topology.select("resid " + str(q_res))
                             if traj_ob.topology.atom(a_idx).element.symbol != 'H' and traj_ob.topology.atom(a_idx).is_sidechain == True]
        
        for l_idx in xrange(len(lig_names)):
            #Initialize the residue bits for this ligand
            temp_bits = [0]*7
            #Create residue and ligand atom-pairs
            a_pairs_dist = list(product(numpy.array(res_heavy_idx), numpy.array(lig_heavy_indices[l_idx])))
            #Compute all pairwise distances and covert to Angstrom; These distances are in nanometers
            pair_distances = md.compute_distances(traj_ob, a_pairs_dist)*10 #shape (1, no_of_a_pairs_dist)
            #Check If ligand is within Apolar interaction distance
            if numpy.any(pair_distances[pair_distances < apolar_cut]):
                #Check if interaction is apolar
                if IsApolar(traj_ob, numpy.array(a_pairs_dist), pair_distances[0],apolar_cut):
                    #Change the apolar bit
                    temp_bits[0] = 1
            #Check If ligand is within Aromatic interaction distance
            if numpy.any(pair_distances[pair_distances < aro_cut]):
                #Check if interaction is aromatic
                isf2f, ise2f = IsAromatic(traj_ob, lig_top[l_idx], numpy.array(a_pairs_dist)
                                                 , pair_distances[0], aro_cut, min(lig_heavy_indices[l_idx]))
                if isf2f:
                    #Change the f2f bit
                    temp_bits[1] = 1
                if ise2f:
                    #Change the f2f bit
                    temp_bits[2] = 1
            #Check If ligand is within Hbond interaction distance
            if numpy.any(pair_distances[pair_distances < hbond_cut]):
                #Check if interaction is Hbond type
                isprotdon, isprotacc = IsHbond(traj_ob, lig_top[l_idx], numpy.array(a_pairs_dist),
                        pair_distances[0], hbond_cut, min(lig_heavy_indices[l_idx]))
                if isprotdon:
                    #Change the protein donor bit
                    temp_bits[3] = 1
                if isprotacc:
                    #Change the protein acceptor bit
                    temp_bits[4] = 1
            #Check If ligand is within Electrostatic interaction distance
            if numpy.any(pair_distances[pair_distances < elec_cut]):
                #Check if interaction is Electrostatic
                prot_pos, prot_neg = IsElectro(traj_ob, lig_top[l_idx], numpy.array(a_pairs_dist),
                                               pair_distances[0], elec_cut, min(lig_heavy_indices[l_idx]))
                
                if prot_pos:
                    #Change the protein positive bit
                    temp_bits[5] = 1
                if prot_neg:
                    #Change the protein negative bit
                    temp_bits[6] = 1
            
            #append temp_bits to sift7bits
            if not sift7bit.has_key(l_idx):
                sift7bit[l_idx] = temp_bits
            else:
                sift7bit[l_idx].extend(temp_bits)
    #Now Combine the bits for different ligands (also corresponding residue indexes) and return it
    if len(lig_names)> 1:
        combined_sift = [lig_bits for key in sorted(sift7bit) for lig_bits in sift7bit[key]]
        #return combined_sift, list(numpy.array(res_index)+1)*len(sift7bit) #+1 since res_index is 0-based
        return combined_sift, (res_index)*len(sift7bit)
    else:
        return sift7bit[0], res_index        
        
def IsApolar(traj_object, pair_indices, pair_dist, apolar_cut):
    #Get the atom pair indices that are within cutoff
    prospective_pairs = pair_indices[pair_dist < apolar_cut]
    isapolar = False
    #Check if any of the pairs is non_polar_atom pair
    for a_pair in prospective_pairs:
        #Get the atom symbol
        a1_symbol = traj_object.topology.atom(a_pair[0]).element.symbol.lower()
        a2_symbol = traj_object.topology.atom(a_pair[1]).element.symbol.lower()
        if a1_symbol in NON_POLAR_ATOMS and a2_symbol in NON_POLAR_ATOMS:
            isapolar = True
            break
    return isapolar

def IsAromatic(traj_object, lig_object, pair_indices, pair_dist, aro_cut, min_lig_heavy_idx):
    #Get the atom pair indices that are within cutoff
    prospective_pairs = pair_indices[pair_dist < aro_cut]
    ise2f = False #Interaction is edge to face
    isf2f = False #Interaction is face to face
    #Combine all ring indices in ligand 
    all_lig_ring_idx = tuple([r_idx for r_idx in chain(*lig_object.rings)])
    #Check if any of the pairs is aromatic pair
    for a_pair in prospective_pairs:
        #Get the residue name (first index is for protein atom)
        a1_res_name = traj_object.topology.atom(a_pair[0]).residue.name.lower()
        #Check if residue is aromatic type
        if a1_res_name in PROTEIN_AROMATIC_RES:
            #Check if protein atom is in aromatic ring
            a1_atom_name = traj_object.topology.atom(a_pair[0]).name.lower()
            if traj_object.topology.atom(a_pair[0]).is_sidechain and a1_atom_name not in ('cb', 'oh'): # oh for tyr
                #Check if the ligand atom is also aromatic
                if lig_object.property[a_pair[1]-min_lig_heavy_idx]['aromatic']: #a_pair[1]-min_lig_heavy_idx is the corresponding index in bymol object
                    ###Check the angles to see if inteaction is f2f or e2f###
                    #Get the neighbors of the ligand atom
                    a2_lig_neighb = lig_object.neighbors[a_pair[1]-min_lig_heavy_idx]
                    #Get the indices of heavy atom neighbors (directly bonded) of ligand atom that are also in ring
                    a2_lig_ring_neighb = tuple([a_nb for a_nb in a2_lig_neighb if a_nb in all_lig_ring_idx])
                    #Get the coordinates for ligand atom in Angstrom and two of its ring neighbor
                    a2_lig_cord = traj_object.xyz[:,a_pair[1],:][0]*10
                    a2_lig_ring_neighb1_cord = traj_object.xyz[:,a2_lig_ring_neighb[0]+min_lig_heavy_idx,:][0]*10
                    a2_lig_ring_neighb2_cord = traj_object.xyz[:,a2_lig_ring_neighb[1]+min_lig_heavy_idx,:][0]*10
                    #Get the cross product for ligand ring
                    lig_ring_normal = GetCross(numpy.array([a2_lig_cord, a2_lig_ring_neighb1_cord, a2_lig_ring_neighb2_cord]))
                    uni_norm_lig_ring = lig_ring_normal / numpy.sqrt((lig_ring_normal*lig_ring_normal).sum()) #Unit-vector
                    #Get the Neighbor atoms (atom objects) of protein atom
                    a1_atom_neighb = GetNeighbor(traj_object, traj_object.topology.atom(a_pair[0]))
                    #Get the coordinates for protein atom in Angstrom and two of its ring neighbor; Hydrogens don't matter
                    a1_atom_cord = traj_object.xyz[:,a_pair[0],:][0]*10
                    a1_atom_neighb1_cord = traj_object.xyz[:,a1_atom_neighb[0].index,:][0]*10
                    a1_atom_neighb2_cord = traj_object.xyz[:,a1_atom_neighb[1].index,:][0]*10
                    #Get the cross product for protein ring
                    prot_ring_normal = GetCross(numpy.array([a1_atom_cord, a1_atom_neighb1_cord, a1_atom_neighb2_cord]))
                    uni_norm_prot_ring = prot_ring_normal / numpy.sqrt((prot_ring_normal*prot_ring_normal).sum()) #Unit-vector
                    #Get the angle between unit normal vectors
                    normal_thetha_deg = numpy.degrees(numpy.arccos(numpy.dot(uni_norm_lig_ring, uni_norm_prot_ring)))  # In radians
                    
                    if normal_thetha_deg <= 30.0 or normal_thetha_deg >= 150.0:
                        isf2f = True
                    if normal_thetha_deg >= 30.0 and normal_thetha_deg <= 150.0:
                        ise2f = True
                    #If both e2f and f2f have been assigned then no need to check further pairs for this residue
                    if isf2f and ise2f:
                        break
    return isf2f, ise2f

def IsHbond(traj_object, lig_object, pair_indices, pair_dist, hbond_cut, min_lig_heavy_idx):
    #Get the atom pair indices that are within cutoff
    prospective_pairs = pair_indices[pair_dist < hbond_cut]
    hbond_prot_acceptor = False
    hbond_prot_donor = False
    #Check if any of the pairs is polar_atom pair
    for a_pair in prospective_pairs:
        #Get the atom symbol
        a1_symbol = traj_object.topology.atom(a_pair[0]).element.symbol.lower()
        a2_symbol = traj_object.topology.atom(a_pair[1]).element.symbol.lower()
        if a1_symbol in POLAR_ATOMS and a2_symbol in POLAR_ATOMS:
            ###Check if angle criteria is satisfied###
            #Get the Neighbor atoms (atom objects) of protein atom
            a1_atom_neighb = GetNeighbor(traj_object, traj_object.topology.atom(a_pair[0]))
            #Check if protein atom is DONOR
            prot_don_hyd = [a1_neig for a1_neig in a1_atom_neighb if a1_neig.element.symbol == 'H']
            if prot_don_hyd:
                is_prot_donor = True
            else:
                is_prot_donor = False
            #Get the neighbors of the ligand atom
            a2_lig_neighb = lig_object.neighbors[a_pair[1]-min_lig_heavy_idx]
            #Check if ligand atom is DONOR; Returns ligand hydrogen index in traj
            lig_don_hyd = [ a2_neig_idx+min_lig_heavy_idx for a2_neig_idx in a2_lig_neighb if lig_object.property[a2_neig_idx]['symbol'] == 'H']
            if lig_don_hyd:
                is_lig_donor = True
            else:
                is_lig_donor = False
            #Get protein atom coordinate
            a1_atom_cord = traj_object.xyz[:,a_pair[0],:][0]*10
            #Get ligand atom coordinate
            a2_lig_cord = traj_object.xyz[:,a_pair[1],:][0]*10
            #Check angles if protein acceptor and ligand donor
            if not is_prot_donor and is_lig_donor:
                #Go through each neighboring hydrogen atom
                for h_neig in lig_don_hyd:
                    #Get the coordinates
                    lig_hyd_cord = traj_object.xyz[:,h_neig,:][0]*10
                    Thetha_degree = GetAngle(numpy.array([a1_atom_cord, lig_hyd_cord, a2_lig_cord]))
                    if Thetha_degree > 135.0:
                        hbond_prot_acceptor = True
                        break #No need to go through rest of the hydrogens
            #Check angles if protein donor and ligand acceptor
            if is_prot_donor and not is_lig_donor:
                #Go through each neighboring hydrogen atom
                for h_neig in prot_don_hyd:
                    #Get the coordinates
                    prot_hyd_cord = traj_object.xyz[:,h_neig.index,:][0]*10
                    Thetha_degree = GetAngle(numpy.array([a2_lig_cord, prot_hyd_cord, a1_atom_cord]))
                    if Thetha_degree > 135.0:
                        hbond_prot_donor = True
                        break #No need to go through rest of the hydrogens
            #If both donor and acceptor have been assigned then no need to check further pairs for this residue
            if hbond_prot_acceptor and hbond_prot_donor:
                break
    return hbond_prot_donor, hbond_prot_acceptor

def IsElectro(traj_object, lig_object, pair_indices, pair_dist, elec_cut, min_lig_heavy_idx):
    #Get the atom pair indices that are within cutoff
    prospective_pairs = pair_indices[pair_dist < elec_cut]
    prot_pos = False
    prot_neg = False
    #Check if any of the pairs is polar_atom pair
    for a_pair in prospective_pairs:
        #Get the residue name (first index is for protein atom)
        a1_res_name = traj_object.topology.atom(a_pair[0]).residue.name.lower()
        #Check if residue is charged type
        if a1_res_name in PROTEIN_POSITIVE + PROTEIN_NEGATIVE:
            #Get the atom symbol
            a1_symbol = traj_object.topology.atom(a_pair[0]).element.symbol.lower()
            #Check if protein atom is in sidechain and is polar
            if traj_object.topology.atom(a_pair[0]).is_sidechain and a1_symbol in POLAR_ATOMS:
                #Check is ligand atom is polar and has non-zero charge
                a2_symbol = traj_object.topology.atom(a_pair[1]).element.symbol.lower()
                if a2_symbol in POLAR_ATOMS and lig_object.property[a_pair[1]-min_lig_heavy_idx]['charge']:
                    #If protein is negatively charged and ligand is postively charged
                    if a1_res_name in PROTEIN_NEGATIVE and lig_object.property[a_pair[1]-min_lig_heavy_idx]['charge'] > 0:
                        prot_neg = True
                    if a1_res_name in PROTEIN_POSITIVE and lig_object.property[a_pair[1]-min_lig_heavy_idx]['charge'] < 0:
                        prot_pos = True
            if prot_neg and prot_pos:
                break       
    return prot_pos, prot_neg
    
def GetNeighbor(traj_object, atom_name):
    '''
    Returns a tuple of Neighbors of atom (atom object is returned)
    atom_name (in the format ex. TYR86-CE1; mdtraj atom object)
    '''
    each_neighb = [list(each_bond) for each_bond in traj_object.topology.bonds if atom_name in each_bond]
    if not each_neighb:
        raise Exception("Topology Object doesn't contain bond information. Use PDB file to define top in Readtraj")
    neighbs = list(set([all_neighb for all_neighb in chain(*each_neighb)]))
    neighbs.remove(atom_name)
    return tuple(neighbs)


def GetCross(cord_array):
    '''
    Cross product of input 3x3 cord_array
    '''
    vec_12 = cord_array[1] - cord_array[0] #Vector tail at atm1 and head at atm2
    vec_13 = cord_array[2] - cord_array[0] #Vector tail at atm1 and head at atm3
    return numpy.cross(vec_12, vec_13)

def GetAngle(cord_array):
    '''
    Returns the angle formed by atom coordinates in input 3x3 cord_array
    '''
    vec_21 = cord_array[0] - cord_array[1] #Vector tail at atm2 and head at atm1
    vec_31 = cord_array[2] - cord_array[1] #Vector tail at atm2 and head at atm3
    dot = numpy.dot(vec_21, vec_31)
    vec_21mod = numpy.sqrt((vec_21*vec_21).sum()) #length of vector 21
    vec_31mod = numpy.sqrt((vec_31*vec_31).sum())
    return numpy.degrees(numpy.arccos(dot / vec_21mod / vec_31mod))



    
    
    
            
            
       
    
        
            
        
    
    
        
    
    
    
    
    
    
    
