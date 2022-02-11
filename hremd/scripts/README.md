# HREMD simulation instructions

## Software requirements:
- `GROMACS (v2018.8)` patched with `PLUMED (v2.5.4)`. For installation, check here: `https://manual.gromacs.org/documentation/2018.8/download.html` and `https://www.plumed.org/download`.
- AMBER. For installation, check here:` https://ambermd.org/GetAmber.php`.
- Python (v2.7 or 3).
- ParmEd. For installation, check here: `https://github.com/ParmEd/ParmEd`.

## Step 1: System file
- Assumes MD system of a ternary complex (e.g. `SMARCA2:PROTAC:VHL`) is already prepared with AMBER tools (e.g., tleap), and `<filename>.parm7` and `<filename>.rst7` files exist.
- Create a directory `inputs` and copy AMBER system parameter/topology `<filename>.parm7` and `<filename>.rst7` files to `inputs` directory.
- Confirm directories `bin`, `config` and `inputs` with `scripts`, `configuration` and `AMBER` system (`<filename>.parm7` and `<filename>.rst7`) files respectively are present in the working directory.

## Step 2: Run HREMD simulation
- Submit a `SLURM` job submission script: `sbatch bin/00.hremd_run.sh`
### Note: Users need to adjust `SLURM parameters`, `GROMACS/PLUMED`, `AMBER` and `ParmEd executables/modules` according to their high performance computing resources.
- `bin/00.hremd_run.sh` will run prepared MD simulations files for `GROMACS/PLUMED`, and performs `energy minimization`, `NVT/NPT equilibration` and `HREMD simulations`.

## Step 3: Extend HREMD simulation
- If users need to extend HREMD simulation, submit a `SLURM` job: `sbatch bin/01.hremd_run_extend.sh`.
### Note: Users need to adjust `SLURM parameters`, `GROMACS/PLUMED`, `AMBER` and `ParmEd executables/modules` according to their high performance computing resources.
- `bin/01.hremd_run_extend.sh` will extend `HREMD simulation` from last step.

