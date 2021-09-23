"""Distance metrics related to common features of receptor-ligand
processes in molecular systems.

The ReceptorDistance class is an abstract class that provides some
common functionality for normalizing reference states, providing
correct indices of receptor and ligand atoms, and a common image
function.

Subclasses of ReceptorDistance need only implement the
'image_distance' function according to their needs.

The UnbindingDistance is a useful metric for enhancing ligand movement
away from the reference bound state conformation.

The RebindingDistance is a useful metric for enhancing the movement of
a ligand towards a reference state.

"""


import logging

import numpy as np

from wepy.util.util import box_vectors_to_lengths_angles

from geomm.grouping import group_pair
from geomm.superimpose import superimpose
from geomm.rmsd import calc_rmsd
from geomm.centering import center_around

from wepy.resampling.distances.receptor import ReceptorDistance

class RebindingDistance(ReceptorDistance):
    """Distance metric for measuring differences between walker states in
    regards to the RMSDs between ligands.

    Images are produced using the ReceptorDistance.image method. The
    distance between images then is the relative difference between
    the ligand RMSDs to the reference state.
    """


    def __init__(self,
                 ligand_idxs = None,
                 binding_site_idxs = None,
                 native_ligand_idxs = None,
                 native_binding_site_idxs = None,
                 ref_state=None,
                 ):


        ## From the superclass

        # the idxs of the ligand and binding site from the whole state
        self._lig_idxs = ligand_idxs
        self._bs_idxs = binding_site_idxs

        # number of atoms in each
        self._n_lig_atoms = len(self._lig_idxs)
        self._n_bs_atoms = len(self._bs_idxs)

        # the idxs used for the whole image
        self._image_idxs = np.concatenate( (self._lig_idxs, self._bs_idxs) )

        # the idxs of the ligand and binding site within the image
        self._image_lig_idxs = np.arange(self._n_lig_atoms)
        self._image_bs_idxs = np.arange(self._n_lig_atoms, self._n_lig_atoms + self._n_bs_atoms)

        ## New stuff

        # Save native state indices

        self._native_ligand_idxs = (ligand_idxs
                                    if native_ligand_idxs is None
                                    else native_ligand_idxs)

        self._native_bs_idxs = (binding_site_idxs
                                if native_binding_site_idxs is None
                                else native_binding_site_idxs)

        # The image idxs for the native state
        self._native_image_idxs = np.concatenate(
            (
                self._native_ligand_idxs,
                self._native_bs_idxs
            ),
        )

        self._ref_state = ref_state

        self.ref_image = self._ref_state['positions'][self._native_image_idxs]


    def _unaligned_image(self, state, ref_state = False):
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

        # Get correct atomic indicies
        if ref_state:

            lig_idxs = self._native_ligand_idxs
            bs_idxs = self._native_bs_idxs
            image_idxs = self._native_image_idxs


        else:
            lig_idxs = self._lig_idxs
            bs_idxs = self._bs_idxs
            image_idxs = self._image_idxs

        # get the box lengths from the vectors
        box_lengths, box_angles = box_vectors_to_lengths_angles(state['box_vectors'])

        # recenter the protein-ligand complex into the center of the
        # periodic boundary conditions

        # regroup the ligand and protein in together
        grouped_positions = group_pair(state['positions'], box_lengths,
                                       bs_idxs, lig_idxs)

        # then center them around the binding site
        centered_positions = center_around(grouped_positions, bs_idxs)

        # slice these positions to get the image
        state_image = centered_positions[image_idxs]

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

        receptor_image : array of float shape (1,)
            The positions of binding site and ligand after
            preprocessing.

        """

        # get the unaligned image
        state_image = self._unaligned_image(state, ref_state = False)

        # then superimpose it to the reference structure
        sup_image, _, _ = superimpose(self.ref_image, state_image, idxs=self._image_bs_idxs)

        rmsd_image = 1 / calc_rmsd(self.ref_image, sup_image, idxs=self._image_lig_idxs)

        return rmsd_image

    def image_distance(self, image_a, image_b):

        # then we get the absolute value of the reciprocals of these rmsd
        # values
        d = abs(image_a - image_b)

        return d
