#!/usr/bin/env python

from __future__ import division
import mdtraj as md
import numpy
from Bython.Cheminfo import Bymol
from Bython.Structure.configstruc import NON_POLAR_ATOMS, PROTEIN_AROMATIC_RES, POLAR_ATOMS, PROTEIN_POSITIVE, PROTEIN_NEGATIVE, PROTEIN_POLAR_SC
from itertools import product, chain

def SIFT_7bit(traj_ob=None, res_index2=None, res_index=None, aro_cut=4.0, apolar_cut=4.5, hbond_cut=3.5,elec_cut=4.0,verbose=True, use_sc_only=False):
    #Make a list of all protein residue indexes if res_index not defined
    if res_index2 == None:
        res_index2 = [traj_ob.topology.atom(atom_index).residue.index for atom_index in traj_ob.topology.select(' chainid 0 1 and name CA')]

    if res_index == None:
        res_index = [traj_ob.topology.atom(atom_index).residue.index for atom_index in traj_ob.topology.select('chainid 2 3 and name CA')]
           
                  
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
        
        for q2_res in res_index2:
            #Initialize the residue bits for this ligand
            temp_bits = [0]*7

	    #Get the residue heavy atom indices
            if not use_sc_only:
                res2_heavy_idx = [a_idx for a_idx in traj_ob.topology.select("resid " + str(q2_res))
                             if traj_ob.topology.atom(a_idx).element.symbol != 'H']
            else:
                #Do only if residue is not Gly
                if traj_ob.topology.residue(q2_res).name.lower() == 'gly':
                    continue
                res2_heavy_idx = [a_idx for a_idx in traj_ob.topology.select("resid " + str(q2_res))
                             if traj_ob.topology.atom(a_idx).element.symbol != 'H' and traj_ob.topology.atom(a_idx).is_sidechain == True]	

            #Create residue and ligand atom-pairs
            a_pairs_dist = list(product(numpy.array(res_heavy_idx), numpy.array(res2_heavy_idx)))
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
                isf2f, ise2f = IsAromatic(traj_ob, numpy.array(a_pairs_dist)
                                                 , pair_distances[0], aro_cut)
                if isf2f:
                    #Change the f2f bit
                    temp_bits[1] = 1
                if ise2f:
                    #Change the f2f bit
                    temp_bits[2] = 1
            #Check If ligand is within Hbond interaction distance
            if numpy.any(pair_distances[pair_distances < hbond_cut]):
                #Check if interaction is Hbond type
                isprotdon, isprotacc = IsHbond(traj_ob, numpy.array(a_pairs_dist),pair_distances[0], hbond_cut)
                if isprotdon:
                    #Change the protein donor bit
                    temp_bits[3] = 1
                if isprotacc:
                    #Change the protein acceptor bit
                    temp_bits[4] = 1
            #Check If ligand is within Electrostatic interaction distance
            if numpy.any(pair_distances[pair_distances < elec_cut]):
                #Check if interaction is Electrostatic
                prot_pos, prot_neg = IsElectro(traj_ob, numpy.array(a_pairs_dist),
                                               pair_distances[0], elec_cut)
                
                if prot_pos:
                    #Change the protein positive bit
                    temp_bits[5] = 1
                if prot_neg:
                    #Change the protein negative bit
                    temp_bits[6] = 1
            
            #append temp_bits to sift7bits
            if not sift7bit.has_key(q2_res):
                sift7bit[q2_res] = temp_bits
            else:
                sift7bit[q2_res].extend(temp_bits)
    #Now Combine the bits for different ligands (also corresponding residue indexes) and return it
    if len(res_index2)> 1:
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

def IsAromatic(traj_object, pair_indices, pair_dist, aro_cut):
    #Get the atom pair indices that are within cutoff
    prospective_pairs = pair_indices[pair_dist < aro_cut]
    ise2f = False #Interaction is edge to face
    isf2f = False #Interaction is face to face
    
    #Check if any of the pairs is aromatic pair
    for a_pair in prospective_pairs:
        #Get the residue name (first index is for protein atom)
        a1_res_name = traj_object.topology.atom(a_pair[0]).residue.name.lower()
	a2_res_name = traj_object.topology.atom(a_pair[1]).residue.name.lower()
        #Check if residue is aromatic type
        if a1_res_name in PROTEIN_AROMATIC_RES and a2_res_name in PROTEIN_AROMATIC_RES:
            #Check if protein atom is in aromatic ring
            a1_atom_name = traj_object.topology.atom(a_pair[0]).name.lower()
            if traj_object.topology.atom(a_pair[0]).is_sidechain and a1_atom_name not in ('cb', 'oh'): # oh for tyr
                #Check if the ligand atom is also aromatic
                a2_atom_name = traj_object.topology.atom(a_pair[1]).name.lower()
                if traj_object.topology.atom(a_pair[1]).is_sidechain and a2_atom_name not in ('cb', 'oh'): # oh for tyr
                    ###Check the angles to see if inteaction is f2f or e2f###
                    
                    #Get the Neighbor atoms (atom objects) of protein atom
                    a1_atom_neighb = GetNeighbor(traj_object, traj_object.topology.atom(a_pair[0]))
                    #Get the coordinates for protein atom in Angstrom and two of its ring neighbor; Hydrogens don't matter
                    a1_atom_cord = traj_object.xyz[:,a_pair[0],:][0]*10
                    a1_atom_neighb1_cord = traj_object.xyz[:,a1_atom_neighb[0].index,:][0]*10
                    a1_atom_neighb2_cord = traj_object.xyz[:,a1_atom_neighb[1].index,:][0]*10
                    #Get the cross product for protein ring
                    prot_ring_normal = GetCross(numpy.array([a1_atom_cord, a1_atom_neighb1_cord, a1_atom_neighb2_cord]))
                    uni_norm_prot_ring = prot_ring_normal / numpy.sqrt((prot_ring_normal*prot_ring_normal).sum()) #Unit-vector
                    

		    #Get the Neighbor atoms (atom objects) of 2nd protein atom
                    a2_atom_neighb = GetNeighbor(traj_object, traj_object.topology.atom(a_pair[1]))
                    #Get the coordinates for protein atom in Angstrom and two of its ring neighbor; Hydrogens don't matter
                    a2_atom_cord = traj_object.xyz[:,a_pair[1],:][0]*10
                    a2_atom_neighb1_cord = traj_object.xyz[:,a2_atom_neighb[0].index,:][0]*10
                    a2_atom_neighb2_cord = traj_object.xyz[:,a2_atom_neighb[1].index,:][0]*10
                    #Get the cross product for protein ring
                    prot_ring_normal = GetCross(numpy.array([a2_atom_cord, a2_atom_neighb1_cord, a2_atom_neighb2_cord]))
                    uni_norm_prot2_ring = prot_ring_normal / numpy.sqrt((prot_ring_normal*prot_ring_normal).sum()) #Unit-vector


		    #Get the angle between unit normal vectors
                    normal_thetha_deg = numpy.degrees(numpy.arccos(numpy.dot(uni_norm_prot2_ring, uni_norm_prot_ring)))  # In radians
                    
                    if normal_thetha_deg <= 30.0 or normal_thetha_deg >= 150.0:
                        isf2f = True
                    if normal_thetha_deg > 30.0 and normal_thetha_deg < 150.0:
                        ise2f = True
                    #If both e2f and f2f have been assigned then no need to check further pairs for this residue
                    if isf2f and ise2f:
                        break
    return isf2f, ise2f

def IsHbond(traj_object, pair_indices, pair_dist, hbond_cut):
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
            #Check if protein 1 atom is DONOR
            prot_don_hyd = [a1_neig for a1_neig in a1_atom_neighb if a1_neig.element.symbol == 'H']
            if prot_don_hyd:
                is_prot_donor = True
            else:
                is_prot_donor = False
            
            #Get the neighbors of the 2nd protein atom
            #Get the Neighbor atoms (atom objects) of protein atom
            a2_atom_neighb = GetNeighbor(traj_object, traj_object.topology.atom(a_pair[1]))
            #Check if protein 2 atom is DONOR
            prot2_don_hyd = [a2_neig for a2_neig in a2_atom_neighb if a2_neig.element.symbol == 'H']
            if prot2_don_hyd:
                is_prot2_donor = True
            else:
                is_prot2_donor = False
            #Get protein atom coordinate
            a1_atom_cord = traj_object.xyz[:,a_pair[0],:][0]*10
            #Get protein 2 atom coordinate
            a2_atom_cord = traj_object.xyz[:,a_pair[1],:][0]*10
            #Check angles if protein acceptor and protein 2 donor
            if not is_prot_donor and is_prot2_donor:
                #Go through each neighboring hydrogen atom
                for h_neig in prot2_don_hyd:
                    #Get the coordinates
                    prot2_hyd_cord = traj_object.xyz[:,h_neig.index,:][0]*10
                    Thetha_degree = GetAngle(numpy.array([a1_atom_cord, prot2_hyd_cord, a2_atom_cord]))
                    if Thetha_degree > 135.0:
                        hbond_prot_acceptor = True
                        break #No need to go through rest of the hydrogens
            #Check angles if protein donor and protein 2 acceptor
            if is_prot_donor and not is_prot2_donor:
                #Go through each neighboring hydrogen atom
                for h_neig in prot_don_hyd:
                    #Get the coordinates
                    prot_hyd_cord = traj_object.xyz[:,h_neig.index,:][0]*10
                    Thetha_degree = GetAngle(numpy.array([a2_atom_cord, prot_hyd_cord, a1_atom_cord]))
                    if Thetha_degree > 135.0:
                        hbond_prot_donor = True
                        break #No need to go through rest of the hydrogens
            #Check if both ligand and protein atom are donor, then Hbond still possible with one as acceptor and other as donor
            if is_prot_donor and is_prot2_donor:
                #Do Hbond determination for both to see which one is acceptor/donor
                #Assuming protein donor, go through each neighboring hydrogen atom of protein donor
                for h_neig in prot_don_hyd:
                    #Get the coordinates
                    prot_hyd_cord = traj_object.xyz[:,h_neig.index,:][0]*10
                    Thetha_degree = GetAngle(numpy.array([a2_atom_cord, prot_hyd_cord, a1_atom_cord]))
                    if Thetha_degree > 135.0:
                        hbond_prot_donor = True
                        break #No need to go through rest of the hydrogens
                #Assuming ligand donor, go through each neighboring hydrogen atom of ligand donor
                for h_neig in prot2_don_hyd:
                    #Get the coordinates
                    prot2_hyd_cord = traj_object.xyz[:,h_neig.index,:][0]*10
                    Thetha_degree = GetAngle(numpy.array([a1_atom_cord, prot2_hyd_cord, a2_atom_cord]))
                    if Thetha_degree > 135.0:
                        hbond_prot_acceptor = True
                        break #No need to go through rest of the hydrogens
                
            #If both donor and acceptor have been assigned then no need to check further pairs for this residue
            if hbond_prot_acceptor and hbond_prot_donor:
                break
    return hbond_prot_donor, hbond_prot_acceptor

def IsElectro(traj_object, pair_indices, pair_dist, elec_cut):
    #Get the atom pair indices that are within cutoff
    prospective_pairs = pair_indices[pair_dist < elec_cut]
    prot_pos = False
    prot_neg = False
    #Check if any of the pairs is polar_atom pair
    for a_pair in prospective_pairs:
        #Get the residue name (first index is for protein atom)
        a1_res_name = traj_object.topology.atom(a_pair[0]).residue.name.lower()
	a2_res_name = traj_object.topology.atom(a_pair[1]).residue.name.lower()
        #Check if residue is charged type
        if a1_res_name in PROTEIN_POSITIVE + PROTEIN_NEGATIVE and a2_res_name in PROTEIN_POSITIVE + PROTEIN_NEGATIVE:
            #Get the atom symbol
            a1_symbol = traj_object.topology.atom(a_pair[0]).element.symbol.lower()
            #Check if protein atom is in sidechain and is polar
            if traj_object.topology.atom(a_pair[0]).is_sidechain and a1_symbol in POLAR_ATOMS:
                #Check if protein 2 atom is in sidechain and polar
                a2_symbol = traj_object.topology.atom(a_pair[1]).element.symbol.lower()
            	if traj_object.topology.atom(a_pair[1]).is_sidechain and a2_symbol in POLAR_ATOMS:
                    #If protein 1 is negatively charged and protein 2 is postively charged
                    if a1_res_name in PROTEIN_NEGATIVE and a2_res_name in PROTEIN_POSITIVE:
                        prot_neg = True
                    if a1_res_name in PROTEIN_POSITIVE and a2_res_name in PROTEIN_NEGATIVE:
                        prot_pos = True
            if prot_neg and prot_pos:
                break       
    return prot_pos, prot_neg
    
def GetNeighbor(traj_object, atom_name):
    '''
    Returns a tuple of Neighbors of atom (atom object is returned)
    atom_name (in the format ex. TYR86-CE1; mdtraj atom object)
    '''
    #print atom_name
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
    
    
    
            
            
       
    
        
            
        
    
    
        
    
    
    
    
    
    
    
