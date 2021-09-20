
# Degrader Ternary Complex Prediction

## Data

All data is available for download from this [Google Bucket][https://console.cloud.google.com/storage/browser/paperdata]


### Weighted Ensemble Simulations (WES)

The WES data can be found in the subsection: `paperdata/WES_data`

All WES data is stored in [HDF5 format][https://www.hdfgroup.org/solutions/hdf5]. 

The data can be read directly using any HDF5 library and is also
compatible with the [wepy][https://github.com/ADicksonLab/wepy]
library.

#### wepy data

The files in the `paperdata/WES_data` directory that are prefixed by
`wepy` are simulation outputs directly from `wepy`.

The filenames are of the format
`wepy_{run_id}_{replicate_index}.wepy.h5` where the `run_id` values
are zero-padded and can be interpreted as:

| run_id | Description |
|--------+-------------|
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
`westpa` are [WESTPA][https://github.com/westpa/westpa] simulation
outputs converted to a single HDF5 file that is compatible with the
analysis tools of `wepy`.

TODO: rename the IDs to match the sim tracker.

| run_id  | Description |
|---------+-------------|
| 24      |             |
| 26 (25) |             |
| 28 (26) |             |
