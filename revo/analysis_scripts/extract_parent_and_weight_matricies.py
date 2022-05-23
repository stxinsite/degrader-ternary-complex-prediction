import numpy as np
import sys

from wepy.hdf5 import WepyHDF5
from wepy.resampling.decisions.clone_merge import MultiCloneMergeDecision

from wepy.analysis import parents

INPUT_H5_FILE = '6hax_binding_results_REVO.wepy.h5'

job_id = sys.argv[1]

run_idx = 0

wepy_h5 = WepyHDF5(INPUT_H5_FILE, mode = 'r')
wepy_h5.open()

n_walkers = wepy_h5.num_init_walkers(run_idx)
n_cycles = wepy_h5.num_run_cycles(run_idx)

# Make Parent Table
resampling_panel = wepy_h5.run_resampling_panel(run_idx)
parent_panel = parents.parent_panel(MultiCloneMergeDecision, resampling_panel)
parent_matrix = np.array(parents.net_parent_table(parent_panel))

np.savetxt('csv_files/job_{}_parent_matrix.csv'.format(job_id), parent_matrix, delimiter = ',')

weights = []

for walker in range(n_walkers):
    print(walker)
    weights.append(np.array(wepy_h5.h5['runs/0/trajectories/{}/weights'.format(walker)])[:,0])


np.savetxt('csv_files/job_{}_weight_matrix.csv'.format(job_id), np.array(weights).T, delimiter = ',')
