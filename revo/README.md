
We suggest doing this in a conda environment using Python 3.7.

First install requirements from conda.

OpenMM:

```sh
conda install -c conda-forge openmm
```

```sh
conda install -c ambermd pytraj parmed
```

The full listing is given in the `env.yaml` file.

```sh
pip install -r requirements.txt
```

To run the code you will need to have the `stx_wepy` module
importable. Typically by setting the `PYTHONPATH` or running from this
directory, e.g.:

```sh
PYTHONPATH=. python python_scripts/amber_ff_triple_distance_metric.py
```


