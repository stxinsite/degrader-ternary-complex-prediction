#!/bin/bash

# Source GROMACS (GMX) and PLUMED executables.
source /bgfs01/insite/utsab.shrestha/programs/gmx_plumed4/gromacs-2018.8/install_dir/bin/GMXRC.bash
source /bgfs01/insite/utsab.shrestha/programs/gmx_plumed4/plumed-2.5.4/sourceme.sh

# PBC removal/correction.
echo "Protein" "!Water" | gmx_mpi trjconv -s topol0.tpr -f traj_comp0.xtc -n index.ndx -pbc cluster -o cluster-protein.xtc

echo 13 "!Water" | gmx_mpi trjconv -s topol0.tpr -f cluster-protein.xtc -n index.ndx -pbc mol -center -o cluster.mol-center-lig.xtc

echo "!Water" | gmx_mpi trjconv -s topol0.tpr -f cluster.mol-center-lig.xtc -n index.ndx -dump 0 -o 0ps.pdb

rm cluster-protein.xtc *#*
