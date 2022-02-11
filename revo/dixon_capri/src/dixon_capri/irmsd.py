"""
Calculates the IRMSD between two proteins.

This is done specifically for a PROTAC system where we have a target protein
(a protein to be degraded), a ligase protein and a protac binding them together.

The IRMSD is calcuated by:

1) Determining the interface residues from a reference structure.
These are residues are where its heavy atoms are within 10 angstroms
of the other protein. 

2) Align the frame to the reference structure. 
The backbone of the interface residues of the frames are superimposed onto the reference structure.

3) The I-RMSD is the RMSD between the backbone atoms of the interface between a given frame and the
reference structure.


Step 1 will needed to done before calling methods in this file.
"""

import numpy as np

from geomm.rmsd import calc_rmsd

from dixon_capri.align import irmsd_align_frame


def calculate_binary_irmsd(reference_positions,
                           positions,
                           unitcell_lengths,                                                                                                                         
                           target_interface_backbone_idx,                                                                                                   
                           ligase_interface_backbone_idx,
                           interface_backbone_idx):

    """
    This function aligns the structure from the simulation
    to the reference structure and calculates the binary IRMSD.

    Parameters
    ----------

    reference_positions: numpy array shape(n_atoms x 3)
        Atomic positions of the reference structure

    positions: numpy array shape (n_atoms x 3)
        Atomic positions of a given simulation frame

    unitcell_lengths: iterable of len 3
        Lengths for the periodic boundary conditions. Used
        for distance calculations and alignment. 

    target_interface_backbone_idx: iterable
        Atomic indicies for the interface backbone atoms
        of the target protein.

   ligase_interface_backbone_idx: iterable
        Atomic indicies for the interface backbone atoms
        of the ligase protein.

   interface_backbone_idx: iterable
        Atomic indicies for the interface backbone atoms
        of both the target and the ligase proteins

    Outputs
    _______

    irmsd: float
        The rmsd between the interface backbone between the simulation 
        frame and the reference structure.

    """

    # Align frames
    aligned_interface_backbone_pos = irmsd_align_frame(reference_positions,
                                                       positions,                                                                                         
                                                       unitcell_lengths,                                                                                
                                                       target_interface_backbone_idx,                                                                              
                                                       ligase_interface_backbone_idx,                                                                         
                                                       interface_backbone_idx)                                                                                                                          
    # Calculate the interface RMSD                                                                                                                                                                 
    irmsd = calc_rmsd(reference_positions[interface_backbone_idx],                                                                                             
                      aligned_interface_backbone_pos)

    return irmsd
