#!/bin/bash -x

# Set up simulation environment
source env.sh

# Clean up from previous/ failed runs
rm -rvf *.log traj_segs seg_logs istates west.h5 *.json slurm.* job_logs/ node_logs/ init_files/
mkdir   seg_logs traj_segs istates job_logs node_logs init_files


# Set pointer to bstate and tstate
BSTATE_ARGS="--bstate-file $WEST_SIM_ROOT/bstates/bstates.txt"

# Run w_init
w_init  $BSTATE_ARGS  --segs-per-state 10  --work-manager=threads "$@"  &> init.log

