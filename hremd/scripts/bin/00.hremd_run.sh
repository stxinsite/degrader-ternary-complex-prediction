#!/bin/bash
#SBATCH -N 1
#SBATCH -J hremd_run
#######SBATCH -o hremd_run.out
#######SBATCH -e hremd_run.err
#SBATCH -n 20 #Number of processes in total for the entire job. This should be equal to number of replicas
#SBATCH --ntasks-per-node=20
#SBATCH --gres=gpu:4
#SBATCH -p debug
#SBATCH --qos=restrained
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --no-requeue

module load cuda
module load amber/latest

# Source GROMACS (GMX) and PLUMED executables.
source /bgfs01/insite/utsab.shrestha/programs/gmx_plumed4/gromacs-2018.8/install_dir/bin/GMXRC.bash
source /bgfs01/insite/utsab.shrestha/programs/gmx_plumed4/plumed-2.5.4/sourceme.sh

#===============================================
# VARIABLES BELOW TO BE ADJUSTED/CHANGED BY USER
#===============================================

HremdMaxH=0.2                  #Maximum hours to run HREMD.
                               #Keep "HremdMaxH" always less than the walltime.

# Total number of replicas to be used.
NumReplicas=20

# 1. Scaling the Hamiltonian of all protein residues.					              
AllProteinResidues='$4=="ALA"||$4=="GLY"||$4=="VAL"||$4=="LEU"||$4=="ILE"||$4=="SER"||$4=="THR"||$4=="PRO"||$4=="PHE"||$4=="TYR"||$4=="TRP"||$4=="MET"||$4=="CYS"||$4=="CYM"||$4=="HIE"||$4=="HID"||$4=="HIP"||$4=="GLN"||$4=="ASN"||$4=="GLU"||$4=="ASP"||$4=="LYS"||$4=="ARG"||$4=="ACE"||$4=="NME"'

# 2. Scaling the Hamiltonian of selected residues. E.g., Hamiltonian of 1LEU, 5PRO and 7VAL are scaled below.
SelectedProteinResidues='$3==1&&$4=="LEU"||$3==5&&$4=="PRO"||$3==7&&$4=="VAL"'

# 3. For not tempering any of the protein residue, use
NoProteinResidue=None

# Choose either #1 or #2 or #3 from above to scale the Hamiltonian of protein residues.
TemperedResidues=$AllProteinResidues

# Scaling the Hamiltonian of resname FWZ (ligand or PROTAC) if present.
# If resname of ligand is different, then change accordingly. If ligand is absent, no action needed.
Ligand='$4=="FWZ"'
#===============================================

#=======================================
# VARIABLES BELOW USUALLY NEED NO CHANGE
#=======================================
gro=complex.gro                 #GMX structure file.						      
top=complex.top                 #GMX parameter file.

Tmin=300                        #Effective temperature of lowest rank (unbiased) replica.
				#Room temperature.

Tmax=425                        #Effective temperature of highest rank replica.
                                #Tmax=425K is optimal value based on benchmark on different system.
				#But if needed, feel free to change if needed.

CheckPt=5                       #Write MD checkpoint file every 5 minutes.                            
MinimizeMaxH=2                  #Maximum hours for minimization.                                      
EquilibMaxH=6                 #Maximum hours for equilibration.                                     
RepEx=500                       #Attempt replica exchange every 500 MD steps (1 ps).
#======================================

# Delete if GMX "top", "tpr" and "gro" files exist from previous run.
rm complex.top complex.gro topol* processed_.top

# Copy system parameter/topology (amber) files in the current directory.
cp inputs/*.parm7 complex.parm7
cp inputs/*.rst7 complex.rst7

# Convert amber "parm7" and "rst7" files to GMX "top" and "gro" files.
cp config/amber2gmx.py . 
python amber2gmx.py

# Minimization.
gmx_mpi grompp -f config/mini.mdp -c "$gro" -p "$top" -o em.tpr
mpirun -mca btl sm,self,tcp -np 1 gmx_mpi mdrun -s em.tpr -deffnm em -cpi -cpt "$CheckPt" -maxh "$MinimizeMaxH" 

# Index file for non-Water atoms.
gmx_mpi make_ndx -f em.gro -o index.ndx<<EOF
!"Water"
q
EOF

# Generate non-Water gro file for position restraining.
echo !Water | gmx_mpi trjconv -s em.gro -f em.gro -n index.ndx -dump 0 -o posre.gro

# NVT Equilibration.
gmx_mpi grompp -f config/nvt.mdp -c em.gro -r posre.gro -p "$top" -n index.ndx -o nvt.tpr
mpirun -mca btl sm,self,tcp -np 1 gmx_mpi mdrun -s nvt.tpr -deffnm nvt -cpi -cpt "$CheckPt" -maxh "$EquilibMaxH"

# NPT Equilibration.
gmx_mpi grompp -f config/npt.mdp -c nvt.gro -r posre.gro -t nvt.cpt -p "$top" -n index.ndx -o npt.tpr
mpirun -mca btl sm,self,tcp -np 1 gmx_mpi mdrun -s npt.tpr -deffnm npt -cpi -cpt "$CheckPt" -maxh "$EquilibMaxH"

# Standard MD tpr file.
gmx_mpi grompp -f config/md.mdp -c npt.gro -t npt.cpt -p "$top" -n index.ndx -o md.tpr

# Append "_" to hot atomtypes for scaling the Hamiltonian.
awk "{if(\$1!=\";\"&&($TemperedResidues)) print \$1,\$2\"_\",\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11;else if(\$1!=\";\"&&($Ligand)) print \$1,\$2\"_\",\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11;else print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11}" complex.top > processed_.top

# Make replicas with geometric progression.
list=$(
awk -v n=$NumReplicas \
    -v tmin=$Tmin \
    -v tmax=$Tmax \
  'BEGIN{for(i=0;i<n;i++){
    t=tmin*exp(i*log(tmax/tmin)/(n-1));
    printf(t); if(i<n-1)printf(",");
  }
}'
)

for ((i=0;i<NumReplicas;i++))
do
        # Choose lambda as T[0]/T[i]
        # Remember that high temperature is equivalent to low lambda.
        lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i+1))';}')
        # Process topology.
        plumed partial_tempering $lambda < processed_.top > topol"$i".top.tmp1
        # Below "sed" is required if running on NEO.
	# If running on different HPC, remove the line below.
        sed '1,13d' topol"$i".top.tmp1 > topol"$i".top.tmp2
        # Below "sed" are required if using amberff14SB force field in GROMACS when top file is processed with PLUMED.
        sed 's/2C/a2C/g' topol"$i".top.tmp2 > topol"$i".top.tmp3
        sed 's/3C/a3C/g' topol"$i".top.tmp3 > topol"$i".top
        rm topol"$i".top.tmp1 topol"$i".top.tmp2 topol"$i".top.tmp3
 	# Prepare tpr for each replica for HREMD.
	# "-maxwarn" is often needed because box could be charged due to scaling coulomb interaction.
	gmx_mpi grompp -maxwarn 1 -o topol$i.tpr -f config/md.mdp -p topol$i.top -c md.tpr -n index.ndx
done

# Check to make sure that there are enough processes.
if [[ $NumReplicas -gt $SLURM_NPROCS ]]; then
    echo "ERROR: Number of replicas greater than the number of tasks allocated for this job. This will disrupt MPI scheduling."
    exit 1
fi

touch plumed.dat
mpirun -np $NumReplicas -mca btl sm,self,tcp --oversubscribe gmx_mpi mdrun -cpi -cpt $CheckPt -maxh $HremdMaxH -plumed plumed.dat -replex $RepEx -hrex -multi $NumReplicas

rm *#*
