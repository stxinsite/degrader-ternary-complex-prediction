
"""                                                                                                                                                                                                       This script will calculate the IRMSD for every frame in a wepy h5 file.                                                                                                                                   

This is done specifically for a PROTAC system where we have a target protein                                                                                                                              (a protein to be degraded), a ligase protein and a protac binding them together.                                                                                                                          

The IRMSD is calcuated by:                                                                                                                                                                                                                                                                                             
1) Determining the interface residues from a reference structure.                                                                                                                                         These are residues are where its heavy atoms are within 10 angstroms                                                                                                                                      of the other protein.                                                                                                                                                                                                                                                                                                  
2) Align the frame to the reference structure.                                                                                                                                                            The backbone of the interface residues of the frames are superimposed onto the reference structure.                                                                                                                                                                                                                                                                                                    
The output will be a csv file that contains the IRMSD.                                                                                                              
"""                                                                                                                                                                                        
import numpy as np
import mdtraj as mdj                                                                                                             
import sys

from dixon_capri.warhead_rmsd import calculate_warhead_rmsd
from dixon_capri.interface_residues import determine_binary_interface_residues                                                                                                                            
from dixon_capri.irmsd import calculate_binary_irmsd                                                                                                                                                       

                                                                                                             
from wepy.hdf5 import WepyHDF5
from wepy.util.util import traj_box_vectors_to_lengths_angles
                                                                                                             
PROTAC_RESNAME = 'FWZ'
TARGET_PROTEIN_RESID_RANGE = (0,114)
LIGASE_PROTEIN_RESID_RANGE = (115, 263)
IRMSD_DIST_CUTOFF = 1 # nm
                                                                                                             
def strip_hydrogen_from_traj(trajectory):

    top = trajectory.top
                                                                                                             
    h_idx = top.select('element "H"')
                                                                                                             
    no_h_list = []

    for atom in range(trajectory.n_atoms):
                                                                                                             
        if atom not in h_idx:
            no_h_list.append(atom)
                                                                                                             
    no_h_traj = trajectory.atom_slice(no_h_list)
                                                                                                             
    return no_h_traj
                                                                                                             
def calculate_warhead_vector(fields_d, *args, **kwargs):
    print('Hi')
    
    positions = fields_d['positions']                                                                                                                                 
    min_irmsd = fields_d['observables/binary_irmsd']
    box_vectors = fields_d['box_vectors']                                                                                                                             
    unitcell_lengths, unitcell_angles = traj_box_vectors_to_lengths_angles(box_vectors)                                                                               
                                                                                                                                                                      
    n_frames = positions.shape[0]                                                                                                                                     
    n_ref_structures = len(ref_positions_list)                                                                                                                        
                                                                                                                                                                      
    warhead_rmsd_walker = np.zeros(n_frames)                                                                                                                                 
                                                                                                                                                                      
    for frame in range(n_frames):                                                                                                                                                          
        for ref_structure in range(n_ref_structures):                                                                                                                                         
            irmsd = calculate_binary_irmsd(ref_positions_list[ref_structure],                                                                                         
                                           positions[frame][sim_no_h_protein_protac_idx],                                                                             
                                           unitcell_lengths[frame],                                                                                                   
                                           target_interface_backbone_idx_list[ref_structure],                                                                         
                                           ligase_interface_backbone_idx_list[ref_structure],                                                                         
                                           interface_backbone_idx_list[ref_structure])
            if irmsd == min_irmsd[frame]:
                warhead_rmsd_walker[frame] = calculate_warhead_rmsd(ref_positions_list[ref_structure],
                                                                    positions[frame][sim_no_h_protein_protac_idx],
                                                                    unitcell_lengths[frame],
                                                                    target_idx_list[ref_structure],
                                                                    ligase_idx_list[ref_structure],
                                                                    protac_idx_list[ref_structure],
                                                                    target_ca_idx_list[ref_structure],
                                                                    warhead_idx_list[ref_structure])
                break

    print('Done')                                                                                                         
    return warhead_rmsd_walker
                                                                                                             
             
wepy_h5_path = str(sys.argv[1])
simulation_pdb_path = str(sys.argv[2])
ref_structure_pdb_list = str(sys.argv[3])
#irmsd_matrix_path = str(sys.argv[4]) 
                                                                                                                                                   
# Load in simulation pdb
sim_pdb = mdj.load_pdb(simulation_pdb_path)
sim_top = sim_pdb.topology

# Get indicies of proteins and protac of interest for simulation ignoring hydrogens                                                                     
sim_no_h_protein_protac_idx = sim_top.select('protein and not element "H" or resname "FWZ" and not element "H"')                           


# Load in reference structure                                                                                                                                        
ref_structure_path_list = np.loadtxt(ref_structure_pdb_list, dtype = str)                                                                                                    
ref_structure_pdb_list = []                                                                                                                                          
for ref_path_list in ref_structure_path_list:                                                                                                                        
    ref_structure_pdb_list.append(mdj.load_pdb(ref_path_list))                                                                                                       

# Remove hydrogen atoms from all reference structures                                                                                                                
no_ref_h_pdb_list = []                                                                                                                                               
no_ref_h_top_list = []                                                                                                                                               
                                                                                                                                                                     
for ref_structure in ref_structure_pdb_list:                                                                                                                         
    no_ref_h_pdb = strip_hydrogen_from_traj(ref_structure)                                                                                                           
    no_ref_h_top = no_ref_h_pdb.top                                                                                                                                  
    no_ref_h_pdb_list.append(no_ref_h_pdb)                                                                                                                           
    no_ref_h_top_list.append(no_ref_h_top)                                                                                                                           
                                                                                                             

# Get indicies of proteins and protac of interest for reference                                                                                                      
target_idx_list = []
target_ca_idx_list = []
ligase_idx_list = []                                                                                                                                                 
protac_idx_list = []                                                                                                                                                 
protein_protac_idx_list = []
warhead_idx_list = []
                                                                                                                                                                     
for no_h_ref_top in no_ref_h_top_list:                                                                                                                               
    
    target_idx = no_h_ref_top.select('resid {} to {}'.format(TARGET_PROTEIN_RESID_RANGE[0], TARGET_PROTEIN_RESID_RANGE[1]))
    target_ca_idx = no_h_ref_top.select('resid {} to {} and name "CA"'.format(TARGET_PROTEIN_RESID_RANGE[0], TARGET_PROTEIN_RESID_RANGE[1]))
    ligase_idx = no_h_ref_top.select('resid {} to {}'.format(LIGASE_PROTEIN_RESID_RANGE[0], LIGASE_PROTEIN_RESID_RANGE[1]))                                          
    protac_idx = no_h_ref_top.select('resname "{}"'.format(PROTAC_RESNAME))
    warhead_idx = no_h_ref_top.select('name "C30" or name "C31" or name "C32" or name "C33" or name "C34" or name "C35" or name "C36" or name "C37" or name "C38" or name "C39" or name "C40" or name "C41" or name "C42" or name "C43" orname "N5" or name "N6" or name "N7" or name "N8" or name "N9" or name "O6"')          
    protein_protac_idx = np.concatenate((target_idx, ligase_idx, protac_idx))                                                                                        
    target_idx_list.append(target_idx)
    target_ca_idx_list.append(target_ca_idx)
    ligase_idx_list.append(ligase_idx)                                                                                                                               
    protac_idx_list.append(protac_idx)                                                                                                                               
    protein_protac_idx_list.append(protein_protac_idx)
    warhead_idx_list.append(warhead_idx)
    
target_resids = np.arange(TARGET_PROTEIN_RESID_RANGE[0], TARGET_PROTEIN_RESID_RANGE[1])                                                                              
ligase_resids = np.arange(LIGASE_PROTEIN_RESID_RANGE[0], LIGASE_PROTEIN_RESID_RANGE[1])                                                                              

# Remove waters from pdb                                                                                                                                              
no_water_ref_pdb_list = []                                                                                                                                            
no_water_ref_top_list = []                                                                                                                                            
                                                                                                                                                                      
for ref_structure in range(len(no_ref_h_pdb_list)):                                                                                                                   
    no_water_ref_pdb = no_ref_h_pdb_list[ref_structure].atom_slice(protein_protac_idx_list[ref_structure])                                                            
    no_water_ref_top = no_water_ref_pdb.topology                                                                                                                      
                                                                                                                                                                      
    no_water_ref_pdb_list.append(no_water_ref_pdb)                                                                                                                    
    no_water_ref_top_list.append(no_water_ref_top)                                                                                                                    
                                                                                                                                                                      
                                                                                                                                                                      
                                                                                                                                                                      
target_interface_resid_list = []                                                                                                                                      
ligase_interface_resid_list = []                                                                                                                                      
                                                                                                                                                                      
for ref_structure in range(len(no_water_ref_pdb_list)):                                                                                            
    target_interface_resid, ligase_interface_resid = determine_binary_interface_residues(no_water_ref_pdb_list[ref_structure],                                       
                                                                                         target_resids,                                                              
                                                                                         target_idx_list[ref_structure],                                             
                                                                                         ligase_resids,                                                              
                                                                                         ligase_idx_list[ref_structure],                                             
                                                                                         IRMSD_DIST_CUTOFF)                                                          
    
    target_interface_resid_list.append(target_interface_resid)                                                                                                       
    ligase_interface_resid_list.append(ligase_interface_resid)                                                                                                       


# Get backbone atomic indicies for interface residues. They will be what                                                                                                                                  
# we use for alignment and calculating the IRMSD from.                                                                                                                                                    

target_interface_backbone_idx_list = []                                                                                                                                                                   
ligase_interface_backbone_idx_list = []                                                                                                                                                                   
interface_backbone_idx_list = []                                                                                                                                                                         
                                                                                                                                                                                                          
for structure in range(len(no_ref_h_pdb_list)):                                                                                                                                                           
    for resid in range(len(target_interface_resid)):                                                                                                                                                      
        if resid == 0:                                                                                                                                                                                    
            target_interface_backbone_idx = no_water_ref_top.select('resid {} and backbone'.format(target_interface_resid[resid]))                                                                        
            
        else:                                                                                                                                                                                             
            target_interface_backbone_idx = np.concatenate((target_interface_backbone_idx, no_water_ref_top.select('resid {} and backbone'.format(target_interface_resid[resid]))))                    
    for resid in range(len(ligase_interface_resid)):                                                                                                                                                      
        if resid == 0:
            ligase_interface_backbone_idx = no_water_ref_top.select('resid {} and backbone'.format(ligase_interface_resid[resid]))                                                           
        else:                                                                                                                                                                                             
            ligase_interface_backbone_idx = np.concatenate((ligase_interface_backbone_idx, no_water_ref_top.select('resid {} and backbone'.format(ligase_interface_resid[resid]))))                       
            
    target_interface_backbone_idx_list.append(target_interface_backbone_idx)                                                                                                                              
    ligase_interface_backbone_idx_list.append(ligase_interface_backbone_idx)                                                                                                                              
    interface_backbone_idx_list.append(np.concatenate((target_interface_backbone_idx, ligase_interface_backbone_idx)))                                                                                    
    
ref_positions_list = []                                                                                                                                                                                   

for ref_structure in no_water_ref_pdb_list:                                                                                                                                                               
    ref_positions_list.append(ref_structure.xyz[0])                                                                                                                                                       


# Load in hdf5 file                                                                                                             
wepy_h5 = WepyHDF5(wepy_h5_path, 'r+')

print('Starting warhead calculations')
# Calculate irmsd
with wepy_h5:    
    warhead_rmsd_matrix = wepy_h5.compute_observable(calculate_warhead_vector,
                                                     ['positions',
                                                      'box_vectors',
                                                      'observables/binary_irmsd'],
                                                     (ref_positions_list,
                                                      target_interface_backbone_idx_list,
                                                      ligase_interface_backbone_idx_list,
                                                      interface_backbone_idx_list,
                                                      target_idx_list,
                                                      ligase_idx_list,
                                                      protac_idx_list,
                                                      target_ca_idx_list,
                                                      warhead_idx_list,
                                                      sim_no_h_protein_protac_idx),
                                                     save_to_hdf5='warhead_rmsd',
                                                     return_results=True)
                                                                                                             
# Save IRMSD Matrix
#np.savetxt(warhead_matrix_path, warhead_rmsd_matrix, delimiter = ',')
