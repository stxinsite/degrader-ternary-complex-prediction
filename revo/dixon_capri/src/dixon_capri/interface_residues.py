"""
Functions used to determine the interface residues.
"""

import mdtraj as mdj
import numpy as np

def calculate_residue_minimum_distance(trajectory,
                                       prot_1_idx,
                                       prot_2_idx,
                                       h_idx):
                                       
    """
    Determines the minimum distance between sets of heavy atoms.
    If hydrogens are included in the atomic indices list they will
    be removed.
    
    Parameters
    ----------
    
    trajectory: MDTRAJ Trajectory Object.
        Trajectory object that contains the 
        atomic distances of the system of interest.
        
   prot_1_idx: list of ints
       List of atomic indices of interest for the first protein.
       
   prot_2_idx: list of ints
       List of atomic indices of interest for the second protein
    
   h_idx: list of ints
       list of hydrogen atom indices in the system
       
   Outputs
   --------
   
   min_dist: float
       minuimum distance between heavy atoms
       between the atomic indice sets.
    """
                                       
    # Slice out hydrogen atoms
    prot_1_idx_no_h = []
    prot_2_idx_no_h = []
     
    for atom_idx in prot_1_idx:
        if atom_idx not in h_idx:
            prot_1_idx_no_h.append(atom_idx)
            
    for atom_idx in prot_2_idx:
        if atom_idx not in h_idx:
            prot_2_idx_no_h.append(atom_idx)  
                  
    # Generate atom pairs
    atom_pairs = []
    for atom_1 in prot_1_idx_no_h:
        for atom_2 in prot_2_idx_no_h:
            atom_pairs.append((atom_1, atom_2))
           
    # Calculate the atomic distance between sets of atomic indicies
    dists = mdj.compute_distances(trajectory, atom_pairs)
   
    # Determine the minimum distance:
    min_dist = np.min(dists)
   
    return min_dist

def determine_binary_interface_residues(pdb, 
                                        target_resid_list,
                                        target_atom_idx, 
                                        ligase_resid_list,
                                        ligase_atom_idx,
                                        interface_cutoff):

    """
    Determines which residues on the target and ligase that
    are within 10 angstroms of the other protein. Only heavy
    atoms are considered.
    
    Parameters
    ----------
    
    pdb: mdtraj trajectory object.
        The reference structure. Must have
        unitcell lengths and angles.
        
    target_resid_list: iterable of len(n_residues_target) containing ints
        List of residue ids for the target protein
        
    target_atom_resid: iterable of len(n_atoms_target) containing ints
        List of atom indicies for the target protein
        
    ligase_resid_list: iterable of len(n_residues_ligase) containing ints
        List of residue ids for the ligase protein 
 
     ligase_atom_resid: iterable of len(n_atoms_ligase) containing ints
        List of atom indicies for the ligase protein
        
    interface_cutoff: float
        How far each residue can be from the protein to be considered in the 
        interface.
        
    Outputs:
    --------
    
    target_interface_residues_list: list containing ints
        List of residues in the target protein that make up the interface 

    ligase_interface_residues_list: list containing ints
        List of residues in the ligase protein that make up the interface 
    """

    top = pdb.top

    # Get hydrogen atom list so we can exclude them from the distance calculation
    h_idx = top.select('element "H"')
 
    # Find target protein interface   
    target_interface_residue_id = []
    for target_residue in target_resid_list:
        target_residue_atom_idx = top.select('resid {}'.format(target_residue))
        
        min_dist = calculate_residue_minimum_distance(pdb, 
                                                      target_residue_atom_idx, 
                                                      ligase_atom_idx,
                                                      h_idx)
        if min_dist < interface_cutoff:
            target_interface_residue_id.append(target_residue)
         
    # Find ligase protein interface
    ligase_interface_residue_id = []
    for ligase_residue in ligase_resid_list:
        ligase_residue_atom_idx = top.select('resid {}'.format(ligase_residue))
        
        min_dist = calculate_residue_minimum_distance(pdb, 
                                                      ligase_residue_atom_idx, 
                                                      target_atom_idx,
                                                      h_idx)
        if min_dist < interface_cutoff:
            ligase_interface_residue_id.append(ligase_residue)
            
    return target_interface_residue_id, ligase_interface_residue_id

def determine_fnat_residue_contacts(pdb, 
                                    target_resid_list, 
                                    ligase_resid_list,
                                    interface_cutoff):

    """
    Determines which residues are in contact between the target
    and ligase proteins in the reference structure. The residues
    with minimum distances between 5 angstroms. The output
    will be lists of lists of atom pairs. Each entry in the
    outer list is for one contact in the reference structure.
    The inner lists contains atom pairs
    
    
    Parameters
    ----------
    
    pdb: mdtraj trajectory object.
        The reference structure. Must have
        unitcell lengths and angles.
        
    target_resid_list: iterable of ints of len(n_target_residues).
        The residue ids for the target protein
 
    ligase_resid_list: iterable of ints of len(n_ligase_residues).
        The residue ids for the target protein
        
    interface_cutoff: float
        How far each residue can be from the protein to be considered in the 
        interface.
        
    Outputs:
    --------
    
    contact_atom_pairs: list of lists of atom idx pairs
        For each target-protein contact, the all non-hydrogen atom pairs
        between the residues that make up that contact.
    
    
    """

    top = pdb.top

    # Get hydrogen atom list so we can exclude them from the distance calculation
    h_idx = top.select('element "H"')
 
    # Get lists of atom indicies for each residue
    target_residue_idx_list = []
    ligase_residue_idx_list = []
    
    for target_residue in target_resid_list:
        target_residue_idx_list.append(top.select('resid {}'.format(target_residue)))

        
    for ligase_residue in ligase_resid_list:
        ligase_residue_idx_list.append(top.select('resid {}'.format(ligase_residue))) 
 
    # Find target protein interface   
    contact_atom_pairs = []
    for target_residue in range(len(target_residue_idx_list)):
        for ligase_residue in range(len(ligase_residue_idx_list)):

            # Calculate the minimum distance between the residues
            min_dist = calculate_residue_minimum_distance(pdb, 
                                                          target_residue_idx_list[target_residue], 
                                                          ligase_residue_idx_list[ligase_residue],
                                                          h_idx)
                                                          
            # Determine if the target-ligase residue pair forms a contact
            if min_dist < interface_cutoff:
                atom_pairs = []
                for target_atom in range(len(target_residue_idx_list[target_residue])):
                    for ligase_atom in range(len(ligase_residue_idx_list[ligase_residue])):
                        atom_pairs.append((target_residue_idx_list[target_residue][target_atom], 
                                           ligase_residue_idx_list[ligase_residue][ligase_atom]))
                contact_atom_pairs.append(atom_pairs)
                        
    return contact_atom_pairs
