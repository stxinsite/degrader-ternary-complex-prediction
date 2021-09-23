import numpy as np
from collections import defaultdict
import logging
import itertools as it
import time

from geomm.grouping import group_pair
from geomm.superimpose import superimpose
from geomm.rmsd import calc_rmsd
from geomm.centering import center_around

import mdtraj as mdj

from wepy.walker import WalkerState
from wepy.util.util import box_vectors_to_lengths_angles
from wepy.util.mdtraj import json_to_mdtraj_topology

from wepy.boundary_conditions.receptor import ReceptorBC



class RebindingBC(ReceptorBC):
    """Boundary condition for doing re-binding simulations of ligands to a
    receptor.
    Implements the ReceptorBC superclass.
    This boundary condition will warp walkers to a number of initial
    states whenever a walker becomes very close to the native (bound)
    state.
    Thus the choice of the 'initial_states' argument should be walkers
    which are completely unbound (the choice of which are weighted by
    'initial_weight') and the choice of 'native_state' should be of a
    ligand bound to the receptor, e.g. X-ray crystallography or docked
    structure.
    The cutoff for the boundary is an RMSD of the walker to the native
    state which is calculated by first aligning and superimposing the
    entire structure according the atom indices specified in
    'binding_site_idxs', and as the name suggests should correspond to
    some approximation of the binding site of the ligand that occurs
    in the native state. Then the raw RMSD of the native and walker
    ligands is calculated. If this RMSD is less than the 'cutoff_rmsd'
    argument the walker is warped.
    PROGRESS is reported for each walker from this rmsd.
    The BC records are never updated.
    """

    # Records of boundary condition changes (sporadic)
    BC_FIELDS = ReceptorBC.BC_FIELDS + ('native_rmsd_cutoff', )
    """The 'native_rmsd_cutoff' is the cutoff used to determine when
    walkers have re-bound to the receptor, which is defined as the
    RMSD of the ligand to the native ligand bound state, when the
    binding sites are aligned and superimposed.
    """

    BC_SHAPES = ReceptorBC.BC_SHAPES + ((1,), )
    BC_DTYPES = ReceptorBC.BC_DTYPES + (np.float, )

    BC_RECORD_FIELDS = ReceptorBC.BC_RECORD_FIELDS + ('native_rmsd_cutoff', )

    # warping (sporadic)
    WARPING_FIELDS = ReceptorBC.WARPING_FIELDS + ()
    WARPING_SHAPES = ReceptorBC.WARPING_SHAPES + ()
    WARPING_DTYPES = ReceptorBC.WARPING_DTYPES + ()

    WARPING_RECORD_FIELDS = ReceptorBC.WARPING_RECORD_FIELDS + ()

    # progress towards the boundary conditions (continual)
    PROGRESS_FIELDS = ReceptorBC.PROGRESS_FIELDS + ('native_rmsd',)
    PROGRESS_SHAPES = ReceptorBC.PROGRESS_SHAPES + (Ellipsis,)
    PROGRESS_DTYPES = ReceptorBC.PROGRESS_DTYPES + (np.float,)

    PROGRESS_RECORD_FIELDS = ReceptorBC.PROGRESS_RECORD_FIELDS + ('native_rmsd', )
    """Records for the state of this record group.
    The 'native_rmsd' is the is the RMSD of the ligand to the native
    ligand bound state, when the binding sites are aligned and
    superimposed.
    """

    def __init__(self, native_state=None,
                 cutoff_rmsd=None,
                 initial_states=None,
                 initial_weights=None,
                 ligand_idxs=None,
                 binding_site_idxs=None,
                 **kwargs):
        """Constructor for RebindingBC.
        Arguments
        ---------
        native_state : object implementing the State interface
            The reference bound state. Will be automatically centered.
        cutoff_rmsd : float
            The cutoff RMSD for considering a walker bound.
        initial_states : list of objects implementing the State interface
            The list of possible states that warped walkers will assume.
        initial_weights : list of float, optional
            List of normalized probabilities of the initial_states
            provided. If not given, uniform probabilities will be
            used.
        ligand_idxs : arraylike of int
            The indices of the atom positions in the state considered
            the ligand.
        binding_site_idxs : arraylike of int
            The indices of the atom positions in the state considered
            the binding site.
        Raises
        ------
        AssertionError
            If any of the following kwargs are not given:
            native_state, initial_states, ligand_idxs, receptor_idxs.
        """

        super().__init__(initial_states=initial_states,
                         initial_weights=initial_weights,
                         ligand_idxs=ligand_idxs,
                         receptor_idxs=binding_site_idxs,
                         **kwargs)

        # test inputs
        assert native_state is not None, "Must give a native state"
        assert cutoff_rmsd is not None, 'Must give a cutoff rmsd'

        native_state_d = native_state.dict()

        # save the native state and center it around it's binding site
        #native_state_d['positions'] = center_around(native_state['positions'], binding_site_idxs)

        native_state = WalkerState(**native_state_d)

        # save attributes
        self._native_state = native_state
        self._cutoff_rmsd = cutoff_rmsd

    @property
    def native_state(self):
        """The reference bound state to which walkers are compared."""
        return self._native_state

    @property
    def cutoff_rmsd(self):
        """The cutoff RMSD for considering a walker bound."""
        return self._cutoff_rmsd

    @property
    def binding_site_idxs(self):
        """The indices of the atom positions in the state considered the binding site."""

        return self._receptor_idxs

    def _progress(self, walker):
        """Calculate if the walker has bound and provide progress record.
        Parameters
        ----------
        walker : object implementing the Walker interface
        Returns
        -------
        is_bound : bool
           Whether the walker is unbound (warped) or not
        progress_data : dict of str : value
           Dictionary of the progress record group fields
           for this walker alone.
        """

        # first recenter the ligand and the receptor in the walker
        box_lengths, box_angles = box_vectors_to_lengths_angles(walker.state['box_vectors'])
        grouped_walker_pos = group_pair(walker.state['positions'], box_lengths,
                                     self.binding_site_idxs, self.ligand_idxs)

        # center the positions around the center of the binding site
        centered_walker_pos = center_around(grouped_walker_pos, self.binding_site_idxs)

        # superimpose the walker state positions over the native state
        # matching the binding site indices only
        sup_walker_pos, _, _ = superimpose(self.native_state['positions'], centered_walker_pos,
                                 idxs=self.binding_site_idxs)

        # calculate the rmsd of the walker ligand (superimposed
        # according to the binding sites) to the native state ligand
        native_rmsd = calc_rmsd(self.native_state['positions'], sup_walker_pos,
                                idxs=self.ligand_idxs)

        # test to see if the ligand is re-bound
        rebound = False
        if native_rmsd <= self.cutoff_rmsd:
            rebound = True

        progress_data = {'native_rmsd' : native_rmsd}

        return rebound, progress_data


class NewRebindingBC(RebindingBC):

    def __init__(self,
                 native_state=None,
                 cutoff_rmsd=0.2,
                 initial_states=None,
                 initial_weights=None,
                 ligand_idxs=None,
                 binding_site_idxs=None,
                 native_state_ligand_idxs=None,
                 native_state_binding_site_idxs=None,
                 **kwargs):


        super().__init__(native_state=native_state,
                         cutoff_rmsd = cutoff_rmsd,
                         initial_states = initial_states,
                         initial_weights = initial_weights,
                         ligand_idxs = ligand_idxs,
                         binding_site_idxs = binding_site_idxs,
                         **kwargs)

        if type(native_state_ligand_idxs) == type(None):
           self.native_ligand_idxs = ligand_idxs
        else:
            self.native_ligand_idxs = native_state_ligand_idxs

        if type(native_state_binding_site_idxs) == type(None):
            self.native_bs_idxs = binding_site_idxs

        else:
            self.native_bs_idxs = native_state_binding_site_idxs


        # number of atoms of in the ligand and binding site

        self._n_lig_atoms = len(self.ligand_idxs)
        self._n_bs_atoms = len(self.binding_site_idxs)


        # the idxs used for the whole image
        self._image_idxs = np.concatenate( (self.ligand_idxs,
                                            self.binding_site_idxs) )

        self._native_image_idxs = np.concatenate((self.native_ligand_idxs,
                                                  self.native_bs_idxs))

        # the idxs of the ligand and binding site within the image
        self._image_lig_idxs = np.arange(self._n_lig_atoms)

        self._image_bs_idxs = np.arange(self._n_lig_atoms,
                                        self._n_lig_atoms + self._n_bs_atoms)

        self.ref_image = self._native_state['positions'][self._native_image_idxs]



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

        if ref_state:
            bs_idxs = self.native_bs_idxs
            lig_idxs = self.native_ligand_idxs
            image_idxs = self._native_image_idxs

        else:
            bs_idxs = self.binding_site_idxs
            lig_idxs = self.ligand_idxs
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


    def _progress(self, walker):
        """Calculate if the walker has bound and provide progress record.

        Parameters
        ----------
        walker : object implementing the Walker interface

        Returns
        -------
        is_bound : bool
           Whether the walker is unbound (warped) or not

        progress_data : dict of str : value
           Dictionary of the progress record group fields
           for this walker alone.

        """

        # first recenter the ligand and the receptor in the walker
        state_image = self._unaligned_image(walker.state, ref_state=False)

        # superimpose the walker state positions over the native state
        # matching the binding site indices only
        sup_walker_pos, _, _ = superimpose(self.ref_image,
                                           state_image,
                                           idxs=self._image_bs_idxs)

        # calculate the rmsd of the walker ligand (superimposed
        # according to the binding sites) to the native state ligand
        native_rmsd = calc_rmsd(self.ref_image, sup_walker_pos,
                                idxs=self._image_lig_idxs)

        # test to see if the ligand is re-bound
        rebound = False
        if native_rmsd <= self.cutoff_rmsd:
            rebound = True

        progress_data = {'native_rmsd' : native_rmsd}

        return rebound, progress_data
