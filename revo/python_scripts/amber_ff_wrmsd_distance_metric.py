from functools import lru_cache
import itertools as it
import os
import os.path as osp
import pickle
import logging
import multiprocessing as mp
from copy import copy, deepcopy
import numpy as np
import mdtraj as mdj
import pytraj as pt
import parmed as pmd
from scipy.spatial.distance import euclidean
from geomm.centroid import centroid

# OpenMM libraries for setting up simulation objects and loading
# the forcefields
import simtk.openmm.app as omma
import simtk.openmm as omm
import simtk.unit as unit


## Wepy classes

# the simulation manager and work mapper for actually running the simulation
from wepy.sim_manager import Manager
from wepy.work_mapper.task_mapper import TaskMapper


# the runner for running dynamics and making and it's particular
# state class
from wepy.runners.openmm import OpenMMRunner, OpenMMState, OpenMMGPUWorker, OpenMMGPUWalkerTaskProcess, UNIT_NAMES, GET_STATE_KWARG_DEFAULTS
from wepy.walker import Walker

# classes for making the resampler
from stx_wepy.resampling.distances.new_rebinding_distance import RebindingDistance
from stx_wepy.resampling.resamplers.epsilon_variation_loss_revo import EpsilonVariationLossREVOResampler

# A standard Boundary condition object for unbinding
from stx_wepy.boundary_conditions.new_rebinding_bc import NewRebindingBC

# standard reporters
from wepy.reporter.hdf5 import WepyHDF5Reporter
from stx_wepy.reporter.walker_pickle_reporter import WalkersPickleReporter
from stx_wepy.reporter.dashboard import DashboardReporter
from wepy.reporter.openmm import OpenMMRunnerDashboardSection
from wepy.reporter.revo.dashboard import REVODashboardSection
from wepy.reporter.receptor.dashboard import RebindingBCDashboardSection

from wepy.util.mdtraj import json_to_mdtraj_topology, mdtraj_to_json_topology

from wepy.util.util import (
    lengths_and_angles_to_box_vectors,
    box_vectors_to_lengths_angles
)

from wepy.util.json_top import (
    json_top_atom_count,
    json_top_subset,
    json_top_residue_df,
    json_top_atom_df,
)

## PARAMETERS

# control what is retrieved from the OpenMM context, this is slightly
# different from the defaults in wepy which gets everything
GET_STATE_KWARGS = {
    'getPositions' : True,
    'getVelocities' : True,
    'getForces' : True,
    'getEnergy' : True,
    'getParameters' : True,
    'getParameterDerivatives' : False,
    'enforcePeriodicBox' : False,
}


# obvious but I hate magic numbers
N_SPATIAL_DIMS = 3

SPATIAL_UNITS = unit.nanometer


# System Parameters
RIGID_WATER = True
REMOVE_CM_MOTION = True
HYDROGEN_MASS = None



# distance cutoff for non-bonded interactions
NONBONDED_CUTOFF = 1.0 * unit.nanometer
# Monte Carlo Barostat
# pressure to be maintained
PRESSURE = 1.0*unit.atmosphere
# temperature to be maintained
TEMPERATURE = 300.0*unit.kelvin
# frequency at which volume moves are attempted
VOLUME_MOVE_FREQ = 50

# amber file constants
AMBER_RESTART_COORDINATE_UNIT = unit.angstrom
RESTRAINT_PROTAC_LIGASE_FORCE_CONSTANT = 2 * (
    unit.kilocalorie_per_mole / unit.angstrom**2)

# Platform used for OpenMM which uses different hardware computation
# kernels. Options are: Reference, CPU, OpenCL, CUDA.

# CUDA is the best for NVIDIA GPUs
PLATFORM = 'Reference'

# Langevin Integrator
FRICTION_COEFFICIENT = 1/unit.picosecond
# step size of time integrations
STEP_TIME = 0.002*unit.picoseconds

# Distance metric parameters, these are not used in OpenMM and so
# don't need the units associated with them explicitly, so be careful!

# distance from the ligand in the crystal structure used to determine
# the binding site, used to align ligands in the Unbinding distance
# metric
BINDING_SITE_CUTOFF = 0.8 # in nanometers

# the residue id for the ligand so that it's indices can be determined

SELECTION_KEYS = (
    'ligase',
    'target',
    'ternary',
    'protac',
    'warhead',
    'linker',
    'ligand',
    'protac+ligase',
    'all',
)


PROTAC_RESNAME = 'FWZ'

WATER_RESNAME = 'HOH'
CHLORIDE_RESNAME = 'Cl-'
POTASSIUM_RESNAME = 'K+'
ZINC_RESNAME = 'Zn'

SOLVENT_RESNAMES = (
    WATER_RESNAME,
    CHLORIDE_RESNAME,
    POTASSIUM_RESNAME,
    ZINC_RESNAME,
)


LINKER_RESID = 'FWZ'
PROT_RES_ID_RANGE = (0, 115)
LIGASE_RES_ID_RANGE = (116, 264)

PROTAC_WARHEAD_ATOM_NAMES = (
    'C30',
    'C31',
    'C32',
    'C33',
    'C34',
    'C35',
    'C36',
    'C38',
    'C39',
    'C40',
    'C41',
    'C41',
    'C42',
    'C43',
    'N5',
    'N6',
    'N7',
    'N8',
    'N9',
    'O6'
)

PROTAC_LINKER_ATOM_NAMES = (
    'C28',
    'C29',
    'C44',
    'C45',
    'C46',
    'C47',
    'C48',
    'C49'
)

PROTAC_LIGAND_ATOM_NAMES = (
    'C1',
    'C2',
    'C3',
    'C4',
    'C5',
    'C6',
    'C7',
    'C8',
    'C9',
    'C10',
    'C11',
    'C12',
    'C13',
    'C14',
    'C15',
    'C16',
    'C17',
    'C18',
    'C19',
    'C20',
    'C21',
    'C22',
    'C23',
    'C24',
    'C25',
    'C26',
    'C27',
    'O1',
    'O2',
    'O3',
    'O4',
    'O5',
    'N1',
    'N2',
    'N3',
    'N4',
    'F1',
    'F2'
)

# Resampler parameters

# the maximum weight allowed for a walker
PMAX = 0.1
# the minimum weight allowed for a walker
PMIN = 1e-50


# boundary condition parameters

# maximum distance between between any atom of the ligand and any
# other atom of the protein, if the shortest such atom-atom distance
# is larger than this the ligand will be considered unbound and
# restarted in the initial state
CUTOFF_DISTANCE = 1.0 # nm
# reporting parameters

# these are the properties of the states (i.e. from OpenMM) which will
# be saved into the HDF5
SAVE_FIELDS = ('positions', 'box_vectors', 'velocities')
# these are the names of the units which will be stored with each
# field in the HDF5
UNITS = UNIT_NAMES
# this is the frequency to save the full system as an alternate
# representation, the main "positions" field will only have the atoms
# for the protein and ligand which will be determined at run time
ALL_ATOMS_SAVE_FREQ = 10
# we can specify some fields to be only saved at a given frequency in
# the simulation, so here each tuple has the field to be saved
# infrequently (sparsely) and the cycle frequency at which it will be
# saved
SPARSE_FIELDS = (('velocities', 10),
                )

## INPUTS/OUTPUTS

# the inputs directory
inputs_dir = osp.realpath('./amber_inputs')
# the outputs path
outputs_dir = osp.realpath('./outputs')
# make the outputs dir if it doesn't exist
os.makedirs(outputs_dir, exist_ok=True)

# inputs filenames
native_amber_parm7_filename = '6hax.chainAB.complex.wat.parm7'
native_amber_rst7_filename = '6hax.chainAB.complex.wat.rst7'
native_pdb_filename = '6hax.chainAB.complex.pdb'

starting_amber_parm7_filename = 'complex.parm7'
starting_amber_rst7_filename = 'relax-6.rst7'

# outputs
hdf5_filename = '6hax_binding_results_REVO.wepy.h5'
dashboard_filename = '6hax_binding_REVO.dash.txt'
parameter_pickle_filename = 'REVO_parameters.pkl'

# normalize the input paths
native_amber_parm7_path = osp.join(inputs_dir, native_amber_parm7_filename)
native_amber_rst7_path = osp.join(inputs_dir, native_amber_rst7_filename)
native_pdb_path = osp.join(inputs_dir, native_pdb_filename)

starting_amber_parm7_path = osp.join(inputs_dir, starting_amber_parm7_filename)

starting_amber_rst7_path = osp.join(inputs_dir, starting_amber_rst7_filename)


# normalize the output paths
hdf5_path = osp.join(outputs_dir, hdf5_filename)
dashboard_path = osp.join(outputs_dir, dashboard_filename)
parameter_pickle_path = osp.join(outputs_dir, parameter_pickle_filename)

# REVO Parameters
CHAR_DIST = 0.1484330303735257
DIST_EXP = 6
EPSILON = 2 # Nanometers
MERGE_DIST = 0.2 # Nanometer
# the maximum weight allowed for a walker
PMAX = 0.1
# the minimum weight allowed for a walker
PMIN = 1e-50
BC_CUTOFF = 0

# Save parameters into a pickle file.
parameter_dic = {'pmin':PMIN,
                 'pmax':PMAX,
                 'char_dist':CHAR_DIST,
                 'alpha': DIST_EXP,
                 'merge_dist':MERGE_DIST,
                 'epsilon':EPSILON,
                 'bc_coutoff': BC_CUTOFF}

pickle.dump( parameter_dic, open( parameter_pickle_path, "wb" ) )


@lru_cache()
def system_selection_idxs(
        json_top,
        selection_key,
):
    """Given the JSON topology and a key defining a selection for a subset
    of atoms in a reference structure gets the atom indices needed for
    selecting them as well as a properly normalized subset of the
    topology for that selection.

    The order of the indices is considered canonical and is guaranteed
    the same as the JSON topology order.

    Supported selections:

    - ligase
        The ligase protein only.

    - target
        The target protein which is designed against.

    - ternary
        The complex of all three molecules: protac, ligase, target.

    - protac
        The entire PROTAC molecule union of 'linker', 'warhead', and 'ligand' moieties

    - warhead
        Only the part of the PROTAC that is designed for binding to the target.

    - linker
        The linker region of the PROTAC which connects the warhead and the ligand.

    - ligand
        The selection of the PROTAC that is designed for the ligase.

    - protac+ligase
        Union of the 'protac' and 'ligase' selections, in that order.


    Note that the 'warhead', 'linker', and 'ligand' selections are
    mutually exclusive and the union of them is equivalent to the
    'protac' selection.

    Parameters
    ----------

    json_top : str
        JSON string topology for the reference topology.

    selection_key : str
        A key that specifies the selection.

    Return
    ------

    sel_idxs : array of int
        Atom indices for the selection.

    sel_tops : str
        JSON string topology for the selection subset.

    """

    # TODO: fill in the rest of the selections

    if selection_key == 'all':

        sel_idxs = np.array(range(json_top_atom_count(json_top)))

        sel_top = json_top

    elif selection_key == 'protac':

        res_df = json_top_residue_df(json_top)
        atom_df = json_top_atom_df(json_top)

        res_idx = res_df[res_df['name'] == PROTAC_RESNAME]['index'].values[0]

        sel_idxs = atom_df[atom_df['residue_key'] == res_idx]['index'].values

        # make the subset topology
        sel_top = json_top_subset(json_top,
                                  sel_idxs)


    elif selection_key == 'warhead':

        atom_df = json_top_atom_df(json_top)

        sel_idxs = atom_df[atom_df['name'].isin(
            PROTAC_WARHEAD_ATOM_NAMES
        )]['index'].values

        sel_top = json_top_subset(json_top,
                                  sel_idxs)


    elif selection_key == 'ligand':
        atom_df = json_top_atom_df(json_top)

        sel_idxs = atom_df[atom_df['name'].isin(
            PROTAC_LIGAND_ATOM_NAMES
        )]['index'].values

        sel_top = json_top_subset(json_top,
                                  sel_idxs)

    elif selection_key == 'linker':
        atom_df = json_top_atom_df(json_top)

        sel_idxs = atom_df[atom_df['name'].isin(
            PROTAC_LINKER_ATOM_NAMES
        )]['index'].values

        sel_top = json_top_subset(json_top,
                                  sel_idxs)


    elif selection_key == 'target':

        res_df = json_top_residue_df(json_top)
        atom_df = json_top_atom_df(json_top)

        res_id_range = range(*PROT_RES_ID_RANGE)

        # get the IDs of the residues in this section
        res_idxs = res_df[res_df['resSeq'].isin(res_id_range)]['index'].values

        sel_idxs = atom_df[atom_df['residue_key'].isin(res_idxs)]['index'].values

        # make the subset topology
        sel_top = json_top_subset(json_top,
                                  sel_idxs)


    elif selection_key == 'ligase':

        res_df = json_top_residue_df(json_top)
        atom_df = json_top_atom_df(json_top)

        res_id_range = range(*LIGASE_RES_ID_RANGE)

        # get the IDs of the residues in this section
        res_idxs = res_df[res_df['resSeq'].isin(res_id_range)]['index'].values

        sel_idxs = atom_df[atom_df['residue_key'].isin(res_idxs)]['index'].values

        # make the subset topology
        sel_top = json_top_subset(json_top,
                                  sel_idxs)


    elif selection_key == 'ternary':

        protac_idxs, _ = system_selection_idxs(
            json_top,
            'protac',
        )
        ligase_idxs, _ = system_selection_idxs(
            json_top,
            'ligase',
        )
        target_idxs, _ = system_selection_idxs(
            json_top,
            'target',
        )

        sel_idxs = np.array(list(it.chain(protac_idxs, ligase_idxs, target_idxs)))

        sel_top = json_top_subset(json_top,
                                  sel_idx)


    elif selection_key == 'protac+ligase':

        protac_idxs, _ = system_selection_idxs(
            json_top,
            'protac',
        )
        ligase_idxs, _ = system_selection_idxs(
            json_top,
            'ligase',
        )

        sel_idxs = np.array(list(it.chain(protac_idxs, ligase_idxs)))

        sel_top = json_top_subset(json_top,
                                  sel_idxs)

    else:
        raise ValueError(f"Unknown selection key: {selection_key}")


    return sel_idxs, sel_top

def load_amber_files(
        prmtop_path,
        nc_rst_path,
):
    """Load the amber starting files and return a normalized JSON top and
    OpenMM ForceField object.

    Parameters
    ----------

    prmtop_path : str or Path
        Path to the 'prmtop' file, must have a 'parm7' file extension.

    nc_rst_path : str or Path
        Path to the NetCDF 'restart' file with the starting system
        coordinates and box vectors, must have a 'rst7' file
        extension.

    Returns
    -------

    positions : array of shape (N_ATOMS, N_DIMS)
       The positions from the Amber restart file for the molecular
       system.

    velocities : array of shape (N_ATOMS, N_DIMS)
       The initial velocities from the Amber restart file for the
       molecular system.

    box_vectors : array of shape (N_DIMS, N_DIMS)
       The box vectors from the Amber restart file.

    json_top : str
        The topology in JSON format, support for this is in wepy.

    amber_prmtop : AmberPrmtopFile
        Special loader class that can create systems from an amber 'prmtop' file.

    """

    ## Force Field

    prmtop = omma.AmberPrmtopFile(str(prmtop_path))

    ## Topology

    # we use the mdtraj loader as the reference topology
    mdj_top = mdj.load_prmtop(str(prmtop_path))

    # convert the mdj top to a JSON one
    json_top = mdtraj_to_json_topology(mdj_top)

    ## State Values

    # load restart file in parmed
    restart_pmd = pmd.amber.Rst7(str(nc_rst_path))

    # load the restart file in pytraj
    restart_pt = pt.load(
            str(nc_rst_path),
            top=str(prmtop_path),
    )

    # need to make sure that this is converted to out system units

    # get the positions as a single frame array
    positions = restart_pt.xyz[0,:,:] * AMBER_RESTART_COORDINATE_UNIT

    # same for velocities
    #velocities = restart_pt.velocities[0,:,:] * AMBER_RESTART_COORDINATE_UNIT

    # get the unitcells from parmed since pytraj can't read them properly
    box_angles = restart_pmd.box[3:6]
    box_lengths = restart_pmd.box[0:3]

    # convert to box vectors
    box_vectors = lengths_and_angles_to_box_vectors(
        *box_lengths,
        *box_angles,
    )

    # apply the unit
    box_vectors = box_vectors * AMBER_RESTART_COORDINATE_UNIT

    return (
        positions,
        #velocities,
        box_vectors,
        json_top,
        prmtop,
    )

def make_omm_system_components(
        positions,
        #velocities,
        box_vectors,
        amber_prmtop,
        json_top,
):

    # convert positions, velocities, box_vectors to our standard units
    positions = positions.in_units_of(SPATIAL_UNITS)
    #velocities = velocities.in_units_of(SPATIAL_UNITS)
    box_vectors = box_vectors.in_units_of(SPATIAL_UNITS)

    ## OpenMM Components

    # make the integrator
    integrator = omm.LangevinIntegrator(
        TEMPERATURE,
        FRICTION_COEFFICIENT,
        STEP_TIME,
    )


    # make the system with the appropriate parameters
    omm_system = amber_prmtop.createSystem(
        nonbondedMethod=omma.CutoffPeriodic,
        nonbondedCutoff=NONBONDED_CUTOFF,
        constraints=omma.HBonds,
        rigidWater=RIGID_WATER,
        removeCMMotion=REMOVE_CM_MOTION,
        hydrogenMass=HYDROGEN_MASS,
    )

    # add the barostat and other forces to the system
    barostat = omm.MonteCarloBarostat(
        PRESSURE,
        TEMPERATURE,
        VOLUME_MOVE_FREQ,
    )

    omm_system.addForce(barostat)

    # get the selection of the ligand part of the protac
    ligand_idxs, _ = system_selection_idxs(
        json_top,
        'ligand',
    )

    ligase_idxs, _ = system_selection_idxs(
        json_top,
        'ligase',
    )

    equilibrium_distance = compute_com_dist(
        ligand_idxs,
        ligase_idxs,
        positions,
    )

    restraint_force = make_protac_ligase_restraint_force(
        RESTRAINT_PROTAC_LIGASE_FORCE_CONSTANT,
        equilibrium_distance,
        ligand_idxs,
        ligase_idxs,
    )

    omm_system.addForce(restraint_force)


    # make the initial state
    init_state = make_openmm_state(
        omm_system,
        integrator,
        positions,
        box_vectors,
    )

    # make the compatible topology object
    mdj_top = json_to_mdtraj_topology(json_top)
    omm_topology = mdj_top.to_openmm()

    return omm_system, init_state, integrator, omm_topology

def make_openmm_state(
        omm_system,
        integrator,
        positions,
        box_vectors,
):
    """Helper function to generate a wepy compatible State object from
    OpenMM components and positions.

    Parameters
    ----------

    omm_system : openmm System object
        A full parametrized system for OpenMM.

    integrator : openmm Integrator object
        An integrator object, this will be copied and the one passed
        in will not be mutated.

    positions : array of shape (N_ATOMS, N_DIMS)
        Positions that will be saved in the state.

    box_vectors : array of shape (N_DIMS, N_DIMS)
       The box vectors from the Amber restart file.

    Returns
    -------

    init_state : wepy.runners.openmm.OpenMMState object
        A wepy compatible State object wrapping an OpenMM State object.

    """

    # we need to make a context to make a state object, which will be
    # wrapped
    #
    # we explicitly use the reference platform here so that we don't
    # accidentally make a CUDA or OpenCL context, which screws up
    # setting up contexts later in the same interpreter session
    platform = omm.Platform.getPlatformByName('Reference')
    context = omm.Context(
        omm_system,
        copy(integrator),
        platform,
    )

    # set the positions
    context.setPositions(positions)

    # set the box vectors
    context.setPeriodicBoxVectors(*box_vectors)


    # get the data from this context so we have a state to start the
    # simulation with, GET_STATE_KWARG_DEFAULTS are just the defaults
    # from the OpenMM module in wepy, change this if you have strange
    # fields in your statea
    get_state_kwargs = dict(GET_STATE_KWARG_DEFAULTS)
    init_sim_state = context.getState(**get_state_kwargs)

    # make the wepy walker state
    init_state = OpenMMState(init_sim_state)

    return init_state

def first_last_atom(topology, segment_id):

    top_dataframe = topology.to_dataframe()[0]
    segment_id_dataframe = top_dataframe[top_dataframe['segmentID'] == segment_id]

    for index, row in segment_id_dataframe.head(1).iterrows():
        first_atom_idx = row['serial'] - 1

    for index, row in segment_id_dataframe.tail(1).iterrows():
        last_atom_idx = row['serial']

    seg_idx = np.arange(first_atom_idx, last_atom_idx)

    return(seg_idx)

def binding_site_atoms(mdtraj_topology, host_idxs, guest_idxs, coords):

    # make a trajectory to compute the neighbors from
    traj = mdj.Trajectory([coords], mdtraj_topology)

    # selects protein atoms which have less than 8 A from ligand
    # atoms in the crystal structure
    neighbors_idxs = mdj.compute_neighbors(traj, BINDING_SITE_CUTOFF, guest_idxs)

    # selects protein atoms from neighbors list
    binding_selection_idxs = np.intersect1d(neighbors_idxs, host_idxs)

    return binding_selection_idxs

def map_native_idx_to_starting_idx_binding_site_atoms(native_mol_idx, start_mol_idx, native_bs_idx):

    starting_bs_idx = []

    for atom in native_bs_idx:

        for native_atom_idx in native_mol_idx:

            if atom == native_mol_idx[native_atom_idx]:

                starting_bs_idx.append(start_mol_idx[native_atom_idx])

    return(starting_bs_idx)

def make_protac_ligase_restraint_force(
        k_param,
        equilibrium_dist,
        protac_ligand_idxs,
        ligase_idxs,

):

    num_groups = 2

    param_sym = "k"


    expr = f"{param_sym}*(distance(g1,g2)-{equilibrium_dist})^2"

    # expr = f"k*distance(g1,g2)^2"

    force = omm.CustomCentroidBondForce(
        num_groups,
        expr,
    )

    force.addPerBondParameter(param_sym)

    # ALERT: you must explicitly cast the integers of the selection
    # indices to plain old python int values for OpenMM. The error
    # doesn't really make this clear either.
    protac_grp_idx = force.addGroup(
        [int(i) for i in protac_ligand_idxs],
    )

    ligase_grp_idx = force.addGroup(
        [int(i) for i in ligase_idxs],
    )

    # create the specific bond
    force.addBond(
        (protac_grp_idx, ligase_grp_idx),
        (k_param,),
    )

    return force

def compute_com_dist(
        sel_idxs_a,
        sel_idxs_b,
        coords,
):

    com_a = centroid(coords[sel_idxs_a])
    com_b = centroid(coords[sel_idxs_b])

    dist = euclidean(com_a, com_b)

    return dist


def main(n_runs, n_cycles, steps, n_walkers, n_workers=1, seed=None):
    ## Load objects needed for various purposes

    native_positions,  native_box_vectors, native_json_top, native_amber_prmtop = load_amber_files(native_amber_parm7_path,
                                                                                                                     native_amber_rst7_path)



    start_positions, start_box_vectors, start_json_top, start_amber_prmtop = load_amber_files(starting_amber_parm7_path,
                                                                                                                starting_amber_rst7_path)



    start_mdj_top = json_to_mdtraj_topology(start_json_top)


    start_pro_b_idx = start_mdj_top.select('resid 0 to 115')
    start_pro_d_idx = start_mdj_top.select('resid 116  to 264')
    start_linker_idx = start_mdj_top.select('name "C30" or name "C31" or name "C32" or name "C33" or name "C34" or name "C35" or name "C36" or name "C37" or name "C38" or name "C39" or name "C40" or name "C41"  or name "C42" or name "C43" or name "N5" or name "N6" or name "N7" or name "N8" or name "N9" or name "O6"')

    start_guest_idx = start_linker_idx




    # load the crystal structure coordinates
    native_pdb = mdj.load_pdb(native_pdb_path)
    native_top = native_pdb.top

    system, init_state, start_integrator, omm_topology = make_omm_system_components(
        start_positions,
        #start_velocities,
        start_box_vectors,
        start_amber_prmtop,
        start_json_top,
    )

    native_system, native_init_state, native_integrator, native_omm_topology = make_omm_system_components(
        native_positions,
        #native_velocities,
        native_box_vectors,
        native_amber_prmtop,
        native_json_top,
    )


    # set up the OpenMMRunner with the system
    runner = OpenMMRunner(system, start_amber_prmtop, start_integrator, platform='CUDA')



    ## Make the distance Metric
    # get the atoms in the binding site according to the crystal structure

    native_pro_b_idx = native_top.select('resid 0 to 115')
    native_pro_d_idx = native_top.select('resid 116 to 264')
    native_warhead_idx = start_mdj_top.select('name "C30" or name "C31" or name "C32" or name "C33" or name "C34" or name "C35" or name "C36" or name "C37" or name "C38" or name "C39" or name "C40" or name "C41"  or name "C42" or name "C43" or name "N5" or name "N6" or name "N7" or name "N8" or name "N9" or name "O6"')




    native_guest_idx = native_warhead_idx

    native_bs_idx = binding_site_atoms(native_pdb.top,
                                       native_pro_b_idx,
                                       native_guest_idx,
                                       native_pdb.xyz[0])

    start_bs_idx = map_native_idx_to_starting_idx_binding_site_atoms(native_pro_b_idx,
                                                                     start_pro_b_idx,
                                                                     native_bs_idx)


    # make the distance metric with the ligand and binding site
    # indices for selecting atoms for the image and for doing the
    # alignments to only the binding site. All images will be aligned
    # to the reference initial state
    reb_distance = RebindingDistance(ligand_idxs = start_guest_idx,
                                     binding_site_idxs = start_bs_idx,
                                     native_ligand_idxs = native_guest_idx,
                                     native_binding_site_idxs = native_bs_idx,
                                     ref_state = native_init_state)

    ## Make the Boundary Conditions

    # makes ref_traj and selects lingand_atom and protein atom  indices
    # instantiate a revo unbindingboudaryconditiobs
    rbc = NewRebindingBC(native_state=native_init_state,
                         cutoff_rmsd=0.00,
                         initial_states=[init_state for i in range(n_walkers)],
                         ligand_idxs=start_guest_idx,
                         binding_site_idxs=start_bs_idx,
                         native_state_ligand_idxs=native_guest_idx,
                         native_state_binding_site_idxs=native_bs_idx)


    ## Make the resampler

    # make a REVO resampler with default parameters and our
    # distance metric
    resampler = EpsilonVariationLossREVOResampler(distance=reb_distance,
                                                  char_dist = CHAR_DIST,
                                                  init_state=init_state,
                                                  merge_dist=MERGE_DIST,
                                                  dist_exponent= DIST_EXP,
                                                  pmax = PMAX,
                                                  pmin = PMIN,
                                                  bc_condition=rbc,
                                                  epsilon = EPSILON,
                                                  best_prog = 'min')


    ## make the reporters

    pkl_reporter = WalkersPickleReporter(save_dir = outputs_dir,
                                         freq = 1,
                                         num_backups = 2)


    # WepyHDF5
    # make a dictionary of units for adding to the HDF5
    # open it in truncate mode first, then switch after first run
    hdf5_reporter = WepyHDF5Reporter(file_path=hdf5_path, mode='w',
                                     # the fields of the State that will be saved in the HDF5 file
                                     save_fields=SAVE_FIELDS,
                                     # the topology in a JSON format
                                     topology=start_json_top,
                                     # the resampler and boundary
                                     # conditions for getting data
                                     # types and shapes for saving
                                     resampler=resampler,
                                     boundary_conditions=rbc,
                                     # the units to save the fields in
                                     units=dict(UNITS),
                                     # sparse (in time) fields
                                     sparse_fields=dict(SPARSE_FIELDS)
                                     # sparse atoms fields
                                     #main_rep_idxs=np.concatenate((prot_idxs, lig_idxs)),
                                     #all_atoms_rep_freq=ALL_ATOMS_SAVE_FREQ
                                    )
    # Dashboard
    # Runner section
    openmm_dashboard_sec = OpenMMRunnerDashboardSection(runner)

    # Receptor section
    re_bc_dashboard_sec = RebindingBCDashboardSection(rbc)

    # Revo Section
    revo_dashboard_sec = REVODashboardSection(resampler)

    # Combine sections for dashboard reporter
    dashboard_reporter = DashboardReporter(file_path = dashboard_filename,
                                           resampler_dash = revo_dashboard_sec,
                                           runner_dash = openmm_dashboard_sec,
                                           bc_dash = re_bc_dashboard_sec)




    reporters = [pkl_reporter, dashboard_reporter, hdf5_reporter]

    mapper = TaskMapper(walker_task_type=OpenMMGPUWalkerTaskProcess,
                        num_workers=n_workers,
                        platform='CUDA',
                        device_ids=[0,1,2,3,4,5,6,7])


    ## Combine all these parts and setup the simulation manager

    # set up parameters for running the simulation
    # initial weights
    init_weight = 1.0 / n_walkers

    # a list of the initial walkers
    init_walkers = [Walker(init_state, init_weight) for i in range(n_walkers)]

    # Instantiate a simulation manager
    sim_manager = Manager(init_walkers,
                          runner=runner,
                          resampler=resampler,
                          boundary_conditions=rbc,
                          work_mapper=mapper,
                          reporters=reporters)


    ### RUN the simulation
    for run_idx in range(n_runs):
        print("Starting run: {}".format(run_idx))
        sim_manager.run_simulation(n_cycles, steps)
        print("Finished run: {}".format(run_idx))


if __name__ == "__main__":
    import time
    import multiprocessing as mp
    import sys
    import logging

    # needs to call spawn for starting processes due to CUDA not
    # tolerating fork
    # mp.set_start_method('spawn', force=True)
    # mp.log_to_stderr(logging.DEBUG)

    if sys.argv[1] == "--help" or sys.argv[1] == '-h':
        print("arguments: n_runs, n_cycles, n_steps, n_walkers, n_workers")
    else:

        n_runs = int(sys.argv[1])
        n_cycles = int(sys.argv[2])
        n_steps = int(sys.argv[3])
        n_walkers = int(sys.argv[4])
        n_workers = int(sys.argv[5])

        print("Number of steps: {}".format(n_steps))
        print("Number of cycles: {}".format(n_cycles))

        steps = [n_steps for i in range(n_cycles)]

        start = time.time()
        main(n_runs, n_cycles, steps, n_walkers, n_workers)
        end = time.time()

        print("time {}".format(end-start))
