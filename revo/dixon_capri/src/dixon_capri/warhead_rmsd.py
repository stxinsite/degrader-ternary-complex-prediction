"""
Calculates the warhead RMSD.

The warhead is the part of the PROTAC molecule which binds to the target protein.

The warhead-RMSD is calcuated by:

1) Align the frame to the reference structure. 
The structure will be aligned to the target C-alpha atoms on the target protein.

2) The warhead-RMSD is the RMSD between the warhead section of the PROTAC between a given frame and the
reference structure.
"""

import numpy as np

from geomm.rmsd import calc_rmsd

from dixon_capri.align import warhead_rmsd_align_frame


def calculate_warhead_rmsd(reference_positions,
                           positions,
                           unitcell_lengths,                                                                                                                         
                           target_idxs,                                                                                                                   
                           ligase_idxs,
                           protac_idxs,
                           target_ca_idxs,
                           warhead_idxs):

    """
    This function aligns the structure from the simulation
    to the reference structure and calculates the warhead RMSD.

    Parameters
    ----------

    reference_positions: numpy array shape(n_atoms x 3)
        Atomic positions of the reference structure

    positions: numpy array shape (n_atoms x 3)
        Atomic positions of a given simulation frame

    unitcell_lengths: iterable of len 3
        Lengths for the periodic boundary conditions. Used
        for distance calculations and alignment. 

    target_idxs: iterable
        Atomic indicies for the target proteins.

   ligase_idxs: iterable
        Atomic indicies for the ligase protein

   protac_idxs: iterable
        Atomic indicies for the protac.

   target_ca_idxs: iterable
        Atomic indicies for the c alpha atoms on the target protein.

   warhead_idxs: iterable
        Atomic indicies for the warhead atoms of the protac.
   

    Outputs
    _______

    warhead_rmsd: float
        The rmsd between the warhead section of the PROTAC
        of the simulation frame and the reference structure.

    """

    # Align frames
    aligned_positions = warhead_rmsd_align_frame(reference_positions,
                                                 positions,                                                                                                           
                                                 unitcell_lengths,
                                                 target_idxs,
                                                 ligase_idxs,
                                                 protac_idxs,
                                                 target_ca_idxs)                                                                                                                 
    # Calculate the interface RMSD                                                                                                                                                                 
    warhead_rmsd = calc_rmsd(reference_positions[warhead_idxs],                                                                                             
                             aligned_positions[warhead_idxs])

    return warhead_rmsd
