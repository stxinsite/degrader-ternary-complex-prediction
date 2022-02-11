# Ternary Structure Prediction: Rosetta-based protocol

## produce a set of ranked predictions of ternary structures POI:PROTAC:E3

### Required Software

ROSETTA

https://rosettacommons.org/software

RDKIT

http://www.rdkit.org/docs/Install.html#how-to-install-rdkit-with-conda

PYTHON2

e.g.:

```
conda create -n py27 python=2.7
```

PYTHON3

e.g.:

```
conda create -n py3 python=3
```

OPENBABEL

e.g.:

```
conda install -c openbabel openbabel
```

### Problem Statement
Quickly produce reasonable ternary structures of POI:PROTAC:E3 complexes for further analysis or refinement with physics-based simulations.

### Definitions and Notations

POI – Protein Of Interest (in the example case presented below – SMARCA2)

E3 – E3-ligase (in the example case presented below – VHL Elongin B)

Warhead – a conventional part of the PROTAC that has large number of contacts with the POI

Ligand – a conventional part of the PROTAC that has large number of contacts with the E3

Linker – a conventional part of the PROTAC that connects warhead and ligand in the peptide chain

### Assumptions

1. binary structures of E3+ligand and POI+warhead are known;

2. structure/SMILES of a linker and attachment vectors to ligand and warhead are known;

### Scope
Must have:

PDB file with E3+ligand and POI+warhead;

PDB file or a SMILES for the linker;

conda environment with RDKit and Python3;

conda environment with Python2 (to use `molfile_to_params.py` from ROSETTA);

### Not in scope:

experimentation with RosettaDock sampling options;

## Workflow
There are 8 conventional sequential steps in the workflow:

1. preparation of PDB files;

2. pre-packing side-chains in the PDB file that will be used at the docking stage;

3. generation of conformers for linkers and infer of constraints to be used in docking;

4. binary docking;

5. re-ranking predictions from binary docking;

6. assembling the ternary structure;

7. ranking;

9. assessment of predictions;

## Step-by-step instructions
The instructions below are describing an application case for the SMARCA2--VHL complex. The initial poses of SMARCA2+warhead and VHL+ligand are taken from the PDB-file 6HAX

### 0. Set Environment Variables and get the scripts

```
export BABEL=/path/to/bin/obabel

export MOL2PARAMS=/path/to/Rosetta/main/source/scripts/python/public/molfile_to_params.py

export ROSETTA_BIN=/path/to/Rosetta/main/source/bin/

export ROSETTA_DATABASE=/path/to/Rosetta/main/database

export BASE_FOLDER=/path/to/the/base/where/the/scripts/folder/is/located/
```

### 1. Prepare initial PDB files

```
mkdir 001.prep
cd 001.prep
```

get the first biological assembly unit from the PDB 6HAX: name the PDB file 6hax.1.pdb ;

```grep 'FWZ' 6hax.1.pdb > 6hax.1.linker.pdb```

the new file `6hax.1.linker.pdb` will contain all atoms of the PROTAC. By visual inspection of the file `6hax.1.linker.pdb` identify atoms of the warhead and ligand and delete them (i.e., only atoms of the linker should remain in the file); change the chainID for the remaining atoms of the linker from `B` to `Z` ;

continue with the next steps:

```cp 6hax.1.pdb 6hax.1.AX_BY.pdb```

in the file `6hax.1.AX_BY.pdb` delete chains `C` , `D`  and waters;  by comparing files `6hax.1.linker.pdb`  and `6hax.1.AX_BY.pdb` remove the atoms of the linker from `6hax.1.AX_BY.pdb` and change chainIDs of the warhead (from `B` to `X` ) and of the ligand (from `B` to `Y`); change the residue name: for the warhead – from `FWZ` to `1_X`; for the ligand – from `FWZ` to `1_Y`;


```grep '1_X' 6hax.1.AX_BY.pdb > warhead.pdb```

the new file `warhead.pdb` contains only the atoms of the warhead; 

```
obabel -ipdb warhead.pdb -O warhead.sdf
conda activate py27   # activate an environment with Python2.7
python $MOL2PARAMS warhead.sdf -p warhead -n 1_X --clobber --keep-names
```

This will produce two new files: `warhead_0001.pdb` and `warhead.params`.

*NOTE:* Rosetta script `molfile_to_params.py` sometimes reorders names of atoms: the order of atomic entries may differ in both the generated PARAMS-file, the output `<X>_0001.pdb` file and the initial `<X>.pdb` file!  If this happens, then both the `docking_prepack_protocol` and `docking_protocol` will produce ERRORs.  To fix this issue, the following steps must be done:

```
conda deactivate
python ../scripts/fix_order_atoms.py warhead.params warhead_0001.pdb warhead.fixed.params
```

The newly produced file `warhead.fixed.params` will be used further.

Perform the same set of step for the ligand:
```
 grep '1_Y' 6hax.1.AX_BY.pdb > ligand.pdb
```

the new file `ligand.pdb` contains only the atoms of the ligand; 

```
obabel -ipdb ligand.pdb -O ligand.sdf
conda activate py27   ## activate an environment with Python2.7
python $MOL2PARAMS ligand.sdf -p ligand -n 1_Y --clobber --keep-names
```

This will produce two new files: `ligand_0001.pdb` and `ligand.params`.

NOTE: Rosetta script `molfile_to_params.py` sometimes reorders names of atoms: the order of atomic entries may differ in both the generated PARAMS-file, the output `<X>_0001.pdb` file and the initial `<X>.pdb` file!  If this happens, then both the `docking_prepack_protocol` and `docking_protocol` will produce ERRORs.  To fix this issue, the following steps must be done:

```
conda deactivate
python ../scripts/fix_order_atoms.pdb ligand.params ligand_0001.pdb ligand.fixed.params
```

The newly produced file `ligand.fixed.params` will be used further.

in the file `warhead_0001.pdb`
change chainID into `X` 

in the file `ligand_0001.pdb`
change chainID into `Y`

use a text editor to replace atom entries for the warhead and ligand in the file `6hax.1.AX_BY.pdb` by the corresponding entries from the files `warhead_0001.pdb`  and `ligand_0001.pdb` respectively.
The newly edited PDB file `6hax.1.AX_BY.pdb` will be used on the next step of the protocol.


### 2.  Rosetta Docking Prepack Protocol

```
cd ..
mkdir 002.prepack
cd 002.prepack
$ROSETTA_BIN/docking_prepack_protocol.static.linuxgccrelease -database $ROSETTA_DATABASE -s ../001.prep/6hax.1.AX_BY.pdb -use_input_sc -extra_res_fa ../001.prep/warhead.fixed.params ../001.prep/ligand.fixed.params
```

This will produce two files: `6hax.1.AX_BY_0001.pdb` and `score.sc`. The produced PDB file will be used in the binary docking stage.

### 3. Generate Conformations for the linker:

Conformations are generated with the `RDKit`

generate an ensemble of conformations:

```
cd 001.prep 
conda activate openbabel
obabel -ipdb 6hax.1.linker.pdb -O 6hax.1.linker.smi 
cd .. 
mkdir 003.linker 
mkdir 003.linker/0.1A 
conda activate my-rdkit-env
python scripts/generate_conformers.py 001.prep/6hax.1.linker.smi 003.linker/0.1A/ 1000 1000 0.1
```

where the last tree parameters in the line are: 1000 – number of conformers to generate, 1000 – number of attempts to make to generate a conformer, 0.1 – the RMS threshold (in A) that each pair of generate conformers must satisfy (run void python scripts/generate_conformers.py to see the instructions for input parameters).

Analyze the generated conformers: get distances between the specified anchor atoms of the linker.
By visual inspection of the PDB files in `003.linker/0.1A/` identify the two terminal atoms – `N1` and `C1`. Use these names in the following command:

```
python scripts/getDistanceStats.py 003.linker/0.3A/ N1 C1 003.linker/0.3A/linker_dist.stat
```

This will produce the file `003.linker/0.1A/linker_dist.stat` that contains distances between atoms `N1` and  `C1` (in A) for each of the generated conformers. The last line in this file contains the mean , standard deviation and number of conformers – see image below:

[plot](https://github.com/stxinsite/degrader-ternary-complex-prediction/blob/taras-dev/docking/pics/fig-1.png)

Screenshot with the contents of the file `003.linker/0.1A/linker_dist.stat`
To visualize the distribution of distances between end of the linker: edit the file `003.linker/0.3A/linker_dist.stat` – remove the last line and add the following line at the top:

[plot](https://github.com/stxinsite/degrader-ternary-complex-prediction/blob/taras-dev/docking/pics/fig-2.png)


Then, run the following commands:

```
conda activate seaborn
python scripts/seaborn_hist_linker.py 003.linker/0.3A/linker_dist.stat
```
And as a result get the following figure:

[plot](https://github.com/stxinsite/degrader-ternary-complex-prediction/blob/taras-dev/docking/pics/fig-3.png)

The values of mean distance  and standard deviation will be used in the restraints on the next step – binary docking with Rosetta.

### 4. Rosetta Docking Protocol
The types and amount of restraints depend on whether we have any additional information about the interactions between proteins or between a protein and a PROTAC. But we can always set at least one harmonic restraint – between the attachment points of the linker. In the example presented below, this harmonic restraint is represented by the first line in the file `004.dock/constraints.sigmoid.1.txt`.

Set distance restraints for docking:

```
mkdir 004.dock
vim 004.dock/constraints.sigmoid.1.txt
```
add the following lines into the file `004.dock/constraints.sigmoid.1.txt`

```
AtomPair N1 1X C9 1Y HARMONIC 8.9 0.7
SiteConstraint CA 34A B SIGMOID 10.0 2.0
SiteConstraint CA 35A B SIGMOID 10.0 2.0
SiteConstraint CA 36A B SIGMOID 10.0 2.0
SiteConstraint CA 37A B SIGMOID 10.0 2.0
SiteConstraint CA 39A B SIGMOID 10.0 2.0
SiteConstraint CA 40A B SIGMOID 10.0 2.0
SiteConstraint CA 41A B SIGMOID 10.0 2.0
SiteConstraint CA 42A B SIGMOID 10.0 2.0
SiteConstraint CA 43A B SIGMOID 10.0 2.0
SiteConstraint CA 45A B SIGMOID 10.0 2.0
SiteConstraint CA 46A B SIGMOID 10.0 2.0
SiteConstraint CA 47A B SIGMOID 10.0 2.0
SiteConstraint CA 48A B SIGMOID 10.0 2.0
SiteConstraint CA 49A B SIGMOID 10.0 2.0
SiteConstraint CA 81A B SIGMOID 10.0 2.0
SiteConstraint CA 82A B SIGMOID 10.0 2.0
SiteConstraint CA 83A B SIGMOID 10.0 2.0
SiteConstraint CA 84A B SIGMOID 10.0 2.0
SiteConstraint CA 85A B SIGMOID 10.0 2.0
SiteConstraint CA 86A B SIGMOID 10.0 2.0
SiteConstraint CA 87A B SIGMOID 10.0 2.0
SiteConstraint CA 88A B SIGMOID 10.0 2.0
SiteConstraint CA 89A B SIGMOID 10.0 2.0
SiteConstraint CA 90A B SIGMOID 10.0 2.0
SiteConstraint CA 91A B SIGMOID 10.0 2.0
SiteConstraint CA 92A B SIGMOID 10.0 2.0
SiteConstraint CA 93A B SIGMOID 10.0 2.0
SiteConstraint CA 94A B SIGMOID 10.0 2.0
SiteConstraint CA 95A B SIGMOID 10.0 2.0
SiteConstraint CA 116B A SIGMOID 10.0 2.0
SiteConstraint CA 118B A SIGMOID 10.0 2.0
SiteConstraint CA 119B A SIGMOID 10.0 2.0
SiteConstraint CA 120B A SIGMOID 10.0 2.0
SiteConstraint CA 121B A SIGMOID 10.0 2.0
SiteConstraint CA 122B A SIGMOID 10.0 2.0
SiteConstraint CA 123B A SIGMOID 10.0 2.0
SiteConstraint CA 124B A SIGMOID 10.0 2.0
SiteConstraint CA 125B A SIGMOID 10.0 2.0
SiteConstraint CA 126B A SIGMOID 10.0 2.0
SiteConstraint CA 128B A SIGMOID 10.0 2.0
SiteConstraint CA 129B A SIGMOID 10.0 2.0
```

where the line 1  indicates that there will be a harmonic restraint between atom `N1` of the “residue“ 1  of the chain `X`  (i.e., the warhead) and atom `C9` of the “residue“ 1  of the chain `Y`  (i.e., the ligand). The values 8.9 and 0.7 are the mean and standard deviation identified on the previous step.
*NOTE:* the distance between these atoms in the original X-tal structure `6HAX` is 7.9 A – see figure below:


   
The line 2  in file `004.dock/constraints.sigmoid.1.txt`  indicates that there will be a sigmoid-type distance restraint between the CA atom of the residue 34 of the chain `A`  and any heavy atom of the chain `B` . Same logic applies to lines 3--42. These restraints are derived from the residues identified to have the highest change on the deuterium incorporation from the HDX experiment – see figure below (blue color represents the residues with the highest deuterium exchange):                 


The values 10.0  and 2.0  in the `SIGMOID` function are some numbers that were not systematically studied (just a guess that turn out to produce good results).

Binary docking


```
mkdir 004.dock/job_1/
cd 004.dock/job_1/
vim slurm_dock.sh
```

In the file `004.dock/job_1/slurm_dock.sh` insert the following lines:


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name 6HAX_HDX__user.name
#SBATCH --partition=project
#SBATCH --qos=maxjobs
#SBATCH --output=constr.1.job1.out
#SBATCH --error=constr.1.job1.err
 
NUM_STR="1000"
$ROSETTA_BIN/docking_protocol.static.linuxgccrelease -database $ROSETTA_DATABASE -s $BASE_FOLDER/002.prepack/6hax.1.AX_BY_0001.pdb -nstruct $NUM_STR -in:file:extra_res_fa  $BASE_FOLDER/001.prep/warhead.fixed.params  $BASE_FOLDER/001.prep/ligand.fixed.params -use_input_sc -docking -dock_pert 2.7 15 -partners AX_BY -ex1 -ex2aro -constraints:cst_file  $BASE_FOLDER/004.dock/constraints.sigmoid.1.txt -constraints:cst_fa_weight 10 -out:file:scorefile score.sc > cst1.job1.mpi.log
wait
echo
echo "All done."
```

Then, submit to SLURM:
`sbatch  slurm_dock.sh` this will produce `$NUM_STR` of docking predictions in the PDB files saved into the current folder (i.e., into `004.dock/job_1/`)

### 5. Re-ranking of Rosetta binary docking predictions 
Use the generated `score.sc` file from the previous stage to re-rank predictions based on the Rosetta interface score, `I_sc`:

```
$ python scripts/getTopNModelsFromScore.py 
please provide:
(i). input score.sc file;
(ii). output score.sc file to save results;
(iii). output file for list of PDB files from Rosetta;
(iv). path to folder where the PDB files are located;

$ python scripts/getTopNModelsFromScore.py 004.dock/job_1/score.sc 004.dock/job_1/score.reranked.sc 004.dock/job_1/listDecoysPDBFiles.txt $BASE_FOLDER/004.dock/job_1/
```

This will produce a new file `score.reranked.sc` with the following line entries separated by commas (i.e., columns will be separated by the ,):

```
# total_score   rms     CAPRI_rank      Fnat    I_sc    I_rms   filename
10.474          9.888   3.000           0.727   -0.025  1.842   6hax.1.AX_BY_0001_0907
-498.495        14.111  1.000           0.273   -0.027  4.023   6hax.1.AX_BY_0001_0501
-486.560        8.286   1.000           0.182   -0.031  5.716   6hax.1.AX_BY_0001_0263
-503.652        10.639  2.000           0.545   -0.032  3.525   6hax.1.AX_BY_0001_0663
-501.564        6.046   2.000           0.545   -0.049  1.693   6hax.1.AX_BY_0001_0037
....
```

Note that one may not always be lucky enough to obtain Top-N predictions with high quality: the RosettaDock employs Monte Carlo random walk to generate and accept/discard the relative poses of docking partners. The `CAPRI_rank`, `rms`, `Fnat` and `I_rms`  are calculated with respect to the initial pose that is supplied as input to RosettaDock. Thus, when doing a blind docking experiment, one has to rely on the values of  the `total_score` and `I_sc`  only for selecting the top-N predictions.

The file `004.dock/job_1/listDecoysPDBFiles.txt`  produced at the previous step contains the paths to PDB files with structures that correspond to the above table:

```
/bgfs01/insite/taras.dauzhenka/data/TPD_Rosetta_6HAX_HDX_reproduction/004.dock/job_1_2/6hax.1.AX_BY_0001_0907.pdb
/bgfs01/insite/taras.dauzhenka/data/TPD_Rosetta_6HAX_HDX_reproduction/004.dock/job_1_2/6hax.1.AX_BY_0001_0501.pdb
/bgfs01/insite/taras.dauzhenka/data/TPD_Rosetta_6HAX_HDX_reproduction/004.dock/job_1_2/6hax.1.AX_BY_0001_0263.pdb
/bgfs01/insite/taras.dauzhenka/data/TPD_Rosetta_6HAX_HDX_reproduction/004.dock/job_1_2/6hax.1.AX_BY_0001_0663.pdb
/bgfs01/insite/taras.dauzhenka/data/TPD_Rosetta_6HAX_HDX_reproduction/004.dock/job_1_2/6hax.1.AX_BY_0001_0037.pdb
...
```

Select the top-10 structures: 

```
cp 004.dock/job_1/listDecoysPDBFiles.txt 004.dock/job_1/listDecoysPDBFiles_top10.txt
```

and delete all lines after the 10th.

Prepare the text file with the list of conformers of the linker that will be aligned onto the predictions from binary docking
```
python scripts/getFileListConformers.py
please provide:
(i). folder with PDB files;
(ii). output text file;

python scripts/getFileListConformers.py $BASE_FOLDER/003.linker/0.1A/ $BASE_FOLDER/003.linker/0.1A/listLinkerPDB.txt
```

The resulting output file `listLinkerPDB.txt`  must contain full paths to the PDB files with the linker conformers.

### 6. Assembling the ternary structures: Alignment of the Linker Conformers onto the Docking Models
This section follows the protocol suggested by Nan Bai et al. ( https://pubs.acs.org/doi/10.1021/acs.jcim.0c01451 ) 

6.1. Preparation
Prepare (by visual inspection) the following text files with the lists of atoms from the linker and ligand/warhead that will be used in structural alignment:

Figure – Illustration of the mechanism of matching atoms during the alignment of linker onto the warhead and ligand.
The file `decoy_atom_list.txt`  contains a list of atoms from the warhead (residue ID `1_X`) and ligand (residue ID `1_Y`) that will be used in the alignment of linker.

```
cat decoy_atom_list.txt
1_X N1
1_X C1
1_X C4
1_Y C9
1_Y C10
1_Y C20
```

The file `decoy_atom_list_delete.txt`  contains a list of atoms from the warhead (residue ID `1_X`) and ligand (residue ID `1_Y`) that will be deleted from the output PDB file if alignment succeeds (these atoms will be replaced by the corresponding atoms from the linker – i.e., by the atoms `N1` and `C1` on the figure above):

```
cat decoy_atom_list_delete.txt
1_X N1
1_Y C9
```

The file `linker_atom_list.txt`  contains a list of atoms from a conformer of the linker that will be used in alignment:

```
cat linker_atom_list.txt
UNL N1
UNL C11
UNL C12
UNL C1
UNL C2
UNL C3
```

The file `linker_atom_list_delete.txt`  contains a list of atoms from the linker that will be deleted from the output PDB file if alignment succeeds (these atoms will be replaced by the corresponding atoms of the warhead (`C1` and `C4` in the Figure above) and the ligand (`C10` and `C20` in the Figure above)):

```
cat linker_atom_list_delete.txt
UNL C11
UNL C12
UNL C2
UNL C3
```
Once these four text files are prepared, the alignments can be run for the whole set of conformers of the linker.

6.2. Alignment

```
mkdir 005.ternary
cp decoy_atom_list.txt 005.ternary
cp decoy_atom_list_delete.txt 005.ternary
cp linker_atom_list.txt 005.ternary
cp linker_atom_list_delete.txt 005.ternary
mkdir 005.ternary/structures/
cd 005.ternary/structures/
```

```
conda activate rdkit_local

time python $BASE_FOLDER/scripts/PROTAC_ternary-master/ternary_model_prediction.py -da $BASE_FOLDER/005.ternary/decoy_atom_list.txt -la $BASE_FOLDER/005.ternary/linker_atom_list.txt -wd $BASE_FOLDER/005.ternary/decoy_atom_list_delete.txt -ld $BASE_FOLDER/005.ternary/linker_atom_list_delete.txt -dl $BASE_FOLDER/004.dock/job_1/listDecoysPDBFiles_top10.txt -ll $BASE_FOLDER/003.linker/0.1A/listLinkerPDB.txt -c 0.3 -r ./rmsd.cst1.job1_0.1A.txt -t specify > log
```

NOTE: It may happen that there is no alignment of linker conformers onto the top-10 binary docking predictions that satisfy the required RMS-threshold (command-line parameter -c in the above script). In this case, one may either 
(1) get back to Section 3 (Generate Conformations for the Linker) and repeat the procedure with the smaller value for rms-threshold in order to get more conformers; or 
(2) increase the number of top-N selected predictions to 50, 100,...  in the `listDecoysPDBFiles_top10.txt`; or
(3) re-run binary docking (Section 4) with an increased number of decoys to produce (`$NUM_STR`). 

### 7. Ranking
Now, as we have a set of predictions of ternary structures, we need to score them to discard those having clashes and leaving those having favorable contacts. Rosetta requires to generate a new PARAMS file for the whole PROTAC molecule that we now have in the PDB file. Thus, we do:

```
cd 005.ternary/structures
mkdir protac
bash ../../scripts/extractPROTACs.sh $BASE_FOLDER/005.ternary/structures/ $BASE_FOLDER/005.ternary/structures/protac/
```
This will (i) remove Hydrogens; (ii) extract PROTAC molecule from the ternary structure; (iii) convert each PDB file the the PROTAC into a corresponding SDF file. Next,

```
conda activate py27
python ../../scripts/runRosettaMol2ParamsForPROTACs.py protac/
```
This will generate the PARAMS files that correspond to the whole PROTAC molecules.

Next, get a SLURM submission script:

```
$ mkdir $BASE_FOLDER/006.rerank
$ conda deactivate
$ python ../../scripts/getSLURMScore.py 
please provide:
(i). input folder with PDB files;
(ii). input folder with PARAMS files;
(iii). file to save slurm script;
(iv). file name pattern: the first five symbols in the PDB filenames that contain ternary structures;

python ../../scripts/getSLURMScore.py $BASE_FOLDER/005.ternary/structures/protac/ $BASE_FOLDER/005.ternary/structures/protac/ $BASE_FOLDER/006.rerank/slurm_score.sh
```

This will create a file `$BASE_FOLDER/006.rerank/slurm_score.sh`  and also will output an integer number into the terminal (in this particular example – 462). This integer is the number of ternary predictions for which the scoring must be run. Thus, we use this number inside the script to set up the submission array (`vim $BASE_FOLDER/006.rerank/slurm_score.sh`):


Submit the jobs: 

```
$ sbatch --array=[0-461]%10 $BASE_FOLDER/006.rerank/slurm_score.sh
```
