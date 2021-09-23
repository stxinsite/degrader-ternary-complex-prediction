import logging

import numpy as np
import mdtraj as mdj

from wepy.util.util import box_vectors_to_lengths_angles

from geomm.grouping import group_pair
from geomm.superimpose import superimpose
from geomm.centering import center_around


from wepy.resampling.distances.distance import Distance


class ProteinProteinContacts(Distance):
    """Distance metric for measuring differences between walker states in
    regards to the number of protein-protein contacts.

    Images are produced using the ReceptorDistance.image method. The
    distance between images then is the relative difference between
    the number of contacts.

    Strength of the contact between the residues between the two protein
    are determined using a sigmoidal function of form y = -1 / (1 + exp(k*(x - x0)))

    where:
    y is the strength of the residue contact.
    x is the minimum distance between the two residues.
    x0 is the distance cutoff (where we want the strength to be at 50%.
    k is the steepness of the curve.


    Parameters
    ----------

    prot_1_resids : list of int
        A list of the residue IDs that you want to be compared to
        the other protein.

    prot_2_resids : list of int
        A list of the residue IDs on protein 1 that you want to be compared to
        the other protein.

    prot_1_idxs : list of int
        List of atom indices identifying protein 1 in the state.

    prot_2_idxs : list of int
        List of atom indices identifying protein 1 in the state.

    TODO: should just be a topology
    trajectory : mdtraj.Trajectory

    TODO: remove this and the alignment since we don't need it
    native_state : 


    distance_cutoff : float

    k : float

    """


    def __init__(self,
                 prot_1_resids= None,
                 prot_2_resids = None,
                 prot_1_idxs = None,
                 prot_2_idxs = None,
                 trajectory = None,
                 native_state = None,
                 distance_cutoff = 0.5,
                 k = 1.0,
                 **kwargs):

        # Ensure all inputs have been given.
        assert prot_1_resids is not None, 'List of resids in protein 1 must be given.'
        assert prot_2_resids is not None, 'List of resids in protein 2 must be given.'
        assert prot_1_idxs is not None, 'List of interface atomic indicies for protein 1 must be given'
        assert prot_2_idxs is not None, 'List of interface atomic indicies for protein 2 must be given'
        assert trajectory is not None, 'MDJ Trajectory must be given.'
        assert native_state is not None, 'Native state must be given.'

        # The resids for the residues of interest
        self.prot_1_resids = prot_1_resids
        self.prot_2_resids = prot_2_resids
        self.prot_1_idxs = prot_1_idxs
        self.prot_2_idxs = prot_2_idxs
        self.trajectory = trajectory
        self.distance_cutoff = distance_cutoff
        self.k = k

        # number of atoms in each
        self._n_prot_1_idxs = len(self.prot_1_idxs)
        self._n_prot_2_idxs = len(self.prot_2_idxs)

        # the idxs used for the whole image
        self._image_idxs = np.concatenate( (self.prot_1_idxs, self.prot_2_idxs) )

        # Get topology of the image.
        self._image_traj = self.trajectory.atom_slice(self._image_idxs)
        self._image_top = self._image_traj.topology
        self._image_h_idx = self._image_top.select('element H')

        # the idxs of the ligand and binding site within the image
        self._image_prot_1_idxs = np.arange(self._n_prot_1_idxs)
        self._image_prot_2_idxs = np.arange(self._n_prot_1_idxs, self._n_prot_2_idxs + self._n_prot_1_idxs)

        # Get resids for image
        self.image_prot_1_resids = np.arange(len(self.prot_1_resids))
        self.image_prot_2_resids = np.arange(len(self.prot_1_resids), len(self.prot_1_resids) + len(self.prot_2_resids))

        # Get heavy atom indicies for each residue in protein 1
        self.image_prot_1_residue_dic = {}

        for resid in self.image_prot_1_resids:
            no_h_idx = []
            resid_idx = self._image_top.select('resid {}'.format(resid))

            for value in resid_idx:
                if value not in self._image_h_idx:
                    no_h_idx.append(value)

            self.image_prot_1_residue_dic.update({resid:np.array(no_h_idx)})

        # Get heavy atom indicies for each residue in protein 1
        self.image_prot_2_residue_dic = {}

        for resid in self.image_prot_2_resids:
            no_h_idx = []
            resid_idx = self._image_top.select('resid {}'.format(resid))

            for value in resid_idx:
                if value not in self._image_h_idx:
                    no_h_idx.append(value)

            self.image_prot_2_residue_dic.update({resid:np.array(no_h_idx)})


        # The image idxs for the native state
        self.ref_image = self._unaligned_image(native_state)




    def _unaligned_image(self, state):
        """The preprocessing method of states.

        First it groups the binding site and ligand into the same
        periodic box image and then centers the box around their
        mutual center of mass and returns only the positions of the
        binding site and ligand.

        Parameters
        ----------
        state : object implementing WalkerState
            State with 'positions' (Nx3 dims) and 'box_vectors' (3x3
            array) attributes.


        Returns
        -------
        """

        # get the box lengths from the vectors
        box_lengths, box_angles = box_vectors_to_lengths_angles(state['box_vectors'])

        # recenter the protein-ligand complex into the center of the
        # periodic boundary conditions

        # regroup the ligand and protein in together
        grouped_positions = group_pair(state['positions'], box_lengths,
                                       self.prot_1_idxs, self.prot_2_idxs)

        # then center them around the binding site
        centered_positions = center_around(grouped_positions, self.prot_1_idxs)

        # slice these positions to get the image
        state_image = centered_positions[self._image_idxs]

        return state_image


    def image(self, state):
        """Transform a state to a receptor image.
        A receptor image is one in which the binding site and ligand
        are first normalized to be within the same periodic box image
        and at the center of the box, and then only the binding site
        and ligand are retained.

        Parameters
        ----------

        state : object implementing WalkerState
            State with 'positions' (Nx3 dims) and 'box_vectors' (3x3 array)
            attributes.

        Returns
        -------
        receptor_image : array of float
            The positions of binding site and ligand after
            preprocessing.

        """

        # get the unaligned image
        state_image = self._unaligned_image(state)

        # then superimpose it to the reference structure
        sup_image, _, _ = superimpose(self.ref_image, state_image, idxs=self._image_prot_1_idxs)

        # calculate the contact strength
        contact_strength = np.sum(self.determine_contacts(sup_image))

        return contact_strength


    def calculate_contact_strength(self, distance):

        contact_strength = 1 / (1 + np.exp(self.k * (distance - self.distance_cutoff)))

        return contact_strength

    def determine_contacts(self, sup_image):
        """Calculate the contact strength between the associated proteins.

        Parameters
        ----------

        sup_image : array of floats (n_atoms, 3)
            Coordinates of state aligned to reference state.


        Returns
        -------

        contact_strengths : array of float (n_atoms * n_atoms, )
            The strength for each individual contact.

        """

        # OPT: this works but might not be a efficient strategy. Could
        # either use mdtraj.compute_contacts, change strategy for min
        # distance, re-implement in Cython or C.

        # Make a vector for the contact strength
        contact_strengths = np.zeros([len(self.image_prot_1_resids) * len(self.image_prot_2_resids)])

        # Create Trajectory object
        traj = mdj.Trajectory(image, self._image_top)

        # Compute minimum distances between each residue
        counter = 0
        for resid_prot_1 in self.image_prot_1_resids:

            # get the atom indices for this residue
            resid_1_idx = self.image_prot_1_residue_dic[resid_prot_1]

            # same thing for protein 2
            for resid_prot_2 in self.image_prot_2_resids:

                resid_2_idx = self.image_prot_2_residue_dic[resid_prot_2]

                # Get atom pairs to calculate minimum distance
                atom_pair = []
                for prot_1_idx in resid_1_idx:

                    for prot_2_idx in resid_2_idx:

                        atom_pair.append((prot_1_idx, prot_2_idx))

                # Compute the atomic distances
                distances = mdj.compute_distances(traj, atom_pair)
                min_distance = np.min(distances)

                # Compute contact strength
                contact_strengths[counter] = self.calculate_contact_strength(min_distance)

                counter += 1


        return contact_strengths


    def image_distance(self, image_a, image_b):

        # then we get absolute difference of the
        # number of contacts between each image
        d = abs(image_a - image_b)

        return d
