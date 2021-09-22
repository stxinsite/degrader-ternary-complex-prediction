
# Degrader Ternary Complex Prediction

## Data

All data is available for download from this [Google Bucket](https://console.cloud.google.com/storage/browser/paperdata)

### Movies

Movies are in the folder `paperdata/movies`. We recommend using the [VLC player](https://www.videolan.org/vlc/) for the best experience.

### Weighted Ensemble Simulations (WES)

The WES data can be found in the subsection: `paperdata/WES_data`

All WES data is stored in [HDF5 format](https://www.hdfgroup.org/solutions/hdf5). 

The data can be read directly using any HDF5 library and is also
compatible with the [wepy](https://github.com/ADicksonLab/wepy)
library.

#### wepy data

The files in the `paperdata/WES_data` directory that are prefixed by
`wepy` are simulation outputs directly from `wepy`.

The filenames are of the format
`wepy_{run_id}_{replicate_index}.wepy.h5` where the `run_id` values
are zero-padded and can be interpreted as:

| run_id | Description |
|--------|-------------|
|      3 |             |
|      4 |             |
|      5 |             |
|      6 |             |
|      8 |             |
|     10 |             |

And the `replicate_index` is just an index of the identical replicate
simulations run with the same initial conditions.

#### WESTPA Data

The files in the `paperdata/WES_data` directory that are prefixed by
`westpa` are [WESTPA](https://github.com/westpa/westpa) simulation
outputs converted to a single HDF5 file that is compatible with the
analysis tools of `wepy`.

The filenames are of the format `westpa_{run_id}.wepy.h5` where the
`run_id` values are zero-padded and can be interpreted as:

| run_id  | Description |
|---------|-------------|
| 24      |             |
| 26      |             |
| 28      |             |


### Folding@Home (FAH) Simulation Data

There is a single FAH project (project ID: 18018) in the
`paperdata/FAH_data` folder. This data is structured in the same
output format as from the FAH infrastructure.

The structure is like:

- `PROJ{project_id}`
  - `RUN{run_id}`
    - `CLONE{clone_id}`
      - `results{gen_index}`
        - `positions.xtc`
        
Where the `positions.xtc` data contain the actual simulation data in
the GROMACS XTC format which can be read by a number of different
software packages.
