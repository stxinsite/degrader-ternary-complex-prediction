"""
Code to align a simulation frame to a reference structure.

Uses geomm package developed by Alex Dickson and his lab.
Repository: https://github.com/ADicksonLab/geomm

"""

import numpy as np

from geomm.grouping import group_pair
from geomm.superimpose import superimpose
from geomm.centering import center_around

def irmsd_align_frame(ref_coords,
                      coords,
                      unitcell_side_lengths,
                      target_interface_backbone_idx,
                      ligase_interface_backbone_idx,
                      interface_backbone_idx):

    """
    Aligns a given frame to the crystal structure using the backbone atoms of the interface residues.
    
    Parameters
    ----------
    
    ref_coords: numpy array of shape (n_protein_atoms + n_protac_atoms, 3)
       Atomic coordinates from the reference structure.
       
    coords: numpy array of shape (n_protein_atoms + n_protac_atoms, 3)
       Atomic coordinates from a frame of the simulation.
       
    unitcell_side_lengths numpy array of len 3.
       The box side lengths to determine periodicity in the system.   
       
    target_interface_backbone_idx: numpy array of ints.
       List of interface bacbone indicies from the target protein.
       
    ligase_interface_backbone_idx: numpy array of ints.
       List of interface bacbone indicies from the ligase protein.
       
    interface_backbone_idx: numpy array of ints.
       List of interface bacbone indicies. A conbination of the previous two arguments.
       
    Outputs
    -------
    
    sup_image numpy array of shape(n_interface_bacbone_atoms x 3)
        Superimposed atomic positions of the interface backbone atoms.
    """

    grouped_positions = group_pair(coords, 
                                   unitcell_side_lengths,
                                   target_interface_backbone_idx, 
                                   ligase_interface_backbone_idx)

    centered_positions = center_around(grouped_positions, interface_backbone_idx)[interface_backbone_idx]

    sup_image, _, _ = superimpose(ref_coords[interface_backbone_idx], centered_positions)

    return(sup_image)

def warhead_rmsd_align_frame(ref_coords,                                                                          
                             coords,                                                                             
                             unitcell_side_lengths,                                                              
                             target_idx,                                                                         
                             ligase_idx,                                                                         
                             protac_idx,                                                                         
                             target_ca_idx):                                                                     
                                                                                                     
    """                                                                                              
    Aligns a given frame to the crystal structure using the backbone atoms of the interface residues\
.                                                                                                    
                                                                                                     
    Parameters                                                                                       
    ----------                                                                                       
                                                                                                     
    ref_coords: numpy array of shape (n_protein_atoms + n_protac_atoms, 3)                           
       Atomic coordinates from the reference structure.                                              
                                                                                                     
    coords: numpy array of shape (n_protein_atoms + n_protac_atoms, 3)                               
       Atomic coordinates from a frame of the simulation.                                            
                                                                                                     
    unitcell_side_lengths numpy array of len 3.                                                      
       The box side lengths to determine periodicity in the system.                                  
                                                                                                     
    target_idx: numpy array of ints.                                              
       Atom indicies of the target protein.                                   
                                                                                                     
    ligase_idx: numpy array of ints.                                              
       Atom indicies of the ligase protein.                                   
                                                                                                     
    ligase_idx: numpy array of ints.                                              
       Atom indicies of the protac.                                   
       
    target_ca_idx: numpy array of ints.                                              
       Atom indicies of C alpha atoms in the target protein.                                   

    Outputs                                                                                          
    -------                                                                                          
                                                                                                     
    sup_image numpy array of shape(n_interface_bacbone_atoms x 3)
        Superimposed atomic positions of the interface backbone atoms.
    """
    
    grouped_positions = group_pair(coords,                                                           
                                   unitcell_side_lengths,                                            
                                   target_idx,                                                       
                                   ligase_idx)                                                       
                                                                                                     
    grouped_positions_2 = group_pair(grouped_positions,                                              
                                     unitcell_side_lengths,                                          
                                     np.concatenate((target_idx, ligase_idx)),                       
                                     protac_idx)                                                     
                                                                                                     
    centered_positions = center_around(grouped_positions_2, target_ca_idx)                           
                                                                                                     
                                                                                                     
    sup_image, _, _ = superimpose(ref_coords, centered_positions, idxs = target_ca_idx)              
                                                                                                     
    return(sup_image)                                                                                



def binary_fnat_center_frame(coords,                                                                             
                             unitcell_side_lengths,                                                              
                             target_idx,
                             ligase_idx,                                                                                                              
                             target_ca_idx):                                                                     
                                                                                                     
    """                                                                                              
    Aligns a given frame to the crystal structure using the backbone atoms of the interface residues\
.                                                                                                    
                                                                                                     
    Parameters                                                                                       
    ----------                                                                                                                                                               
    coords: numpy array of shape (n_protein_atoms + n_protac_atoms, 3)                               
       Atomic coordinates from a frame of the simulation.                                            
                                                                                                     
    unitcell_side_lengths numpy array of len 3.                                                      
       The box side lengths to determine periodicity in the system.                                  
                                                                                                     
    target_idx: numpy array of ints.                                              
       Atom indicies of the target protein.                                   
                                                                                                     
    ligase_idx: numpy array of ints.                                              
       Atom indicies of the ligase protein.                                   
                                                                                                     
    ligase_idx: numpy array of ints.                                              
       Atom indicies of the protac.                                   
       
    target_ca_idx: numpy array of ints.                                              
       Atom indicies of C alpha atoms in the target protein.                                   
    Outputs                                                                                          
    -------                                                                                          
                                                                                                     
    centered_positions numpy array of shape(n_atoms, 3)
        Positions where the target and ligase are in the same
        periodic cell and are at the center of the cell.
    """
    
    grouped_positions = group_pair(coords,                                                           
                                   unitcell_side_lengths,                                            
                                   target_idx,                                                       
                                   ligase_idx)                                                       
                                                                                                     
                                                                                                     
    centered_positions = center_around(grouped_positions, target_ca_idx)                           
                                                                                                                                                      
    return(centered_positions)                                                                                
