import numpy as np
import mdtraj as mdj

from wepy.hdf5 import WepyHDF5
from wepy.util.util import traj_box_vectors_to_lengths_angles

from wepy.resampling.decisions.clone_merge import MultiCloneMergeDecision
#from wepy.boundary_conditions.rebinding import RebindingBC
from wepy.analysis.parents import parent_panel, net_parent_table, parent_table_discontinuities, ancestors

from geomm.grouping import group_pair
from geomm.superimpose import superimpose
from geomm.rmsd import calc_rmsd
from geomm.centering import center_around


import sys


def first_last_c_a(topology, segment_id):

    top_dataframe = topology.to_dataframe()[0]
    segment_id_dataframe = top_dataframe[top_dataframe['segmentID'] == segment_id]

    for index, row in segment_id_dataframe.head(1).iterrows():
        first_atom_idx = row['serial'] - 1

    for index, row in segment_id_dataframe.tail(1).iterrows():
        last_atom_idx = row['serial']

    seg_idx = np.arange(first_atom_idx, last_atom_idx)

    return(seg_idx)


def element_idx(topology, element):
    """                                                                                                                                                                                                  
    Determines all the atom indicies for a specific element                                                                                                                                                

    Inputs:                                                                                                                                                                                                        topology: mdtraj toplogy                                                                                                                                                                                             topology of system of interest                                                                                                                                                                   element: string                                                                                                                                                                                                     Desired element to determine the indicies of.                                                                                                                                                                                                                                                                                                                                                            Outputs:                                                                                                                                                                                                       element_idx: numpy array                                                                                                                                                                                                array containing the atomic indicies for                                                                                                                                                                   a single element in the protein.                                                                                                                                                                                                                                                                                                                                                                                                                           
    """
    element_idx = []

    protein_idx = topology.select('protein')

    data_frame = topology.to_dataframe()[0]

    for index, row in data_frame.iterrows():
        if row['element'] == element:
            if row['serial'] - 1 in protein_idx:
                element_idx.append(row['serial'] - 1)

    return(np.array(element_idx))


def allign_frame(ref_coords,
                 coords,
                 unitcell_side_lengths,
                 btk_idx,
                 cIAP_idx,
                 btk_ca_idx):


    grouped_positions = group_pair(coords, unitcell_side_lengths,
                                    btk_idx, cIAP_idx)

    centered_positions = center_around(grouped_positions, btk_ca_idx)

    sup_image, _, _ = superimpose(ref_coords, centered_positions, idxs=btk_ca_idx)

    return(sup_image)


# What is the last cycle and walker                  
# to make the trajectory                                                                                                                                                                                                                                       
job_id = sys.argv[1]
cycle_idx = int(sys.argv[2])
walker_idx = int(sys.argv[3])
run_idx = 0



pdb = mdj.load_pdb('../amber_inputs/6hax.chainAB.complex.pdb')
topology = pdb.topology


wepy_h5 = WepyHDF5('copy_6hax_binding_results_REVO.wepy.h5', mode='r')

# Make Parent Matrix                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
with wepy_h5:
    resampling_panel = wepy_h5.run_resampling_panel(run_idx)
    parent_panel = parent_panel(MultiCloneMergeDecision, resampling_panel)
    parent_matrix = net_parent_table(parent_panel)
    #parent_matrix = np.array(parent_table_discontinuities(UnbindingBC, parent_matrix, wepy_h5.warping_records([0])))


# Get mdtraj trajectory for the lineage                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
lineage = ancestors(parent_matrix, cycle_idx, walker_idx)

box_length_list = []
box_angle_list = []
for cycle_walker_id in range(len(lineage)):
    print(cycle_walker_id)
    with wepy_h5:

        traj_data = wepy_h5.get_run_trace_fields(run_idx, [lineage[cycle_walker_id]], ['positions', 'box_vectors'])

        positions = traj_data['positions']

        box_vectors = traj_data['box_vectors']

        unitcell_lengths, unitcell_angles = traj_box_vectors_to_lengths_angles(box_vectors)

        if cycle_walker_id == 0:
            positions_all = np.zeros([len(lineage), positions[0].shape[0], positions[0].shape[1]])

        positions_all[cycle_walker_id] = positions[0]

        box_length_list.append(unitcell_lengths)

        box_angle_list.append(unitcell_angles)


box_length_list = np.array(box_length_list).reshape(len(lineage), 3)
box_angle_list = np.array(box_angle_list).reshape(len(lineage), 3)

#traj = mdj.Trajectory(positions_all,
#                      topology,
#                      unitcell_lengths = box_length_list,
#                      unitcell_angles = box_angle_list)

#traj.save_dcd('job_{}_cycle_{}_walker_{}_trajectory.dcd'.format(job_id, cycle_idx, walker_idx))

#print('Trajectory made')

# Superimpose to BTK c alpha                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
#btk_idx = first_last_c_a(topology, seg_id)
btk_idx = topology.select('resid 0 to 115')
ciap_idx = topology.select('resid 116  to 263')                                                                                                                                             

pro_ca = topology.select('protein and name "CA"')

#ciap_idx = first_last_c_a(topology, 'PROD')
tky_idx = topology.select('resname "FWZ"')

ciap_tky_idx = np.concatenate((ciap_idx, tky_idx))
pro_lig_idx = np.concatenate((btk_idx, ciap_idx, tky_idx))


pro_lig_traj = pdb.atom_slice(pro_lig_idx)
pro_lig_top = pro_lig_traj.top

pro_lig_traj.save_pdb('no_water_SMARCA_2_VHL_bound.pdb')

btk_ca = []
for atom in pro_ca:
    if atom in btk_idx:
        btk_ca.append(atom)

print('Superpose start')

superpose_pos = np.zeros([positions_all.shape[0], len(pro_lig_idx), 3])

for cycle in range(positions_all.shape[0]):
    superpose_pos[cycle] = allign_frame(pdb.xyz[0][pro_lig_idx],
                                        positions_all[cycle][pro_lig_idx],
                                        box_length_list[cycle],
                                        btk_idx,
                                        ciap_tky_idx,
                                        btk_ca)

superpose_traj = mdj.Trajectory(superpose_pos,
                                pro_lig_top,
                                unitcell_lengths = box_length_list,
                                unitcell_angles = box_angle_list)



#superpose_traj = traj.superpose(pdb,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#                                frame = 0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
#                                atom_indices = btk_ca)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               

print('superpose finished')

superpose_traj.save_dcd('job_{}_cycle_{}_walker_{}_SMARC2_superpose_trajectory.dcd'.format(job_id, cycle_idx, walker_idx))
