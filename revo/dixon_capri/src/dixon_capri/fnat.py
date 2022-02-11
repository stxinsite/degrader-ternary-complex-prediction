import mdtraj as mdj
import numpy as np

from dixon_capri.align import binary_fnat_center_frame 


def calculate_binary_fnat(positions,
                          topology,
                          unitcell_lengths,
                          unitcell_angles,
                          contact_atom_pair_list,
                          contact_cutoff):

    """
    Calculates the fraction of native contacts (fnat). A contact is when a residue on the 
    target protein and ligase protein are within a cutoff distance. 5 angstroms is
    typically used as a cutoff. 

    The fnat is calculated by:

    1: Determine which residues make a contact in a reference structure. This step is
    done by the calc_fnat_contacts function found in the interface.py file. 

    2: Ensure the target and ligase protein are within the same periodic cell for your
    simulation frame.

    3: Determine how many of the contacts found in the native structure are present in your
    simulation frame.

    4: fnat = frame contacts / total number of native contacts

    Parameters
    ----------

    positions: numpy array of shape (n_atoms,3)
       Atomic positions of your frame.

    topology: mdtraj topology object
       Topology of the system. Needed to make the mdtraj Trajectory object to calculate distances.

    unitcell_lengths: numpy array of len(3)
       Box lengths of periodic cell.

    unitcell_angles: numpy array of len (3)
       Box angles of periodic cell.

    contact_atom_pair_list: list of lists containing tuple of ints of atom pairs.
       List of each contact made in the reference structure. Each entry contains a
       series of tuples that contain two atom indicies of pairs between all atoms
       of the target protein residue atoms and ligase residue atoms that form a
       contact.

    contact_cutoff: float
       Distance between residues that define if they are in contact.

    Outputs:
    --------

    fnat: float
      Fraction of native contacts formed in a simulation frame

    """
    trajectory = mdj.Trajectory(positions,
                                topology,
                                unitcell_lengths = unitcell_lengths,
                                unitcell_angles = unitcell_angles)

    contacts = 0

    for residue in range(len(contact_atom_pair_list)):

        contacts += determine_contacts(trajectory,
                                       contact_atom_pair_list[residue],
                                       contact_cutoff)

    fnat = contacts / len(contact_atom_pair_list)

    return fnat



def determine_contacts(trajectory,
                       atom_pairs_list,
                       contact_cutoff):

    """
    Determines if the target-ligase residue pair forms a contact.
    Residues will form a contact if th minimum distance of heavy
    atoms is under a cutoff. Generally 5 angstroms 

    Parameters
    ----------

    trajectory: mdtraj Trajectory
        Trajectory of the aligned simulation frame
 
    atom_pairs_list: list of tuples (atom_1, atom_2)
        List of atom pairs betwen heavy atoms on a target residue
        and a ligase residue. These residues form a contact
        in the reference structure.
 
    contact_cutoff: float
        Distance cutoff that determines if two residues
        form a contact. Normally 5 angstroms.
 
    Outputs:
    --------
 
    contact: float
        Float that indicates if the residues form a contact.
        1 indicates a contact is formed.
        0 indicates there is no contact formed.

    """

    distance = mdj.compute_distances(trajectory, atom_pairs_list)

    min_distance = np.min(distance)

    if min_distance < contact_cutoff:
        contact = 1
    else:
        contact = 0

    return contact

                           
