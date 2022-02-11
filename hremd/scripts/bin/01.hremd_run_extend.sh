#!/bin/bash
#SBATCH -N 1
#SBATCH -J hremd_run
######SBATCH -o hremd_run.out
######SBATCH -e hremd_run.err
#SBATCH -n 20 #Number of processes in total for the entire job. This should be equal to number of replicas
#SBATCH --ntasks-per-node=20 
#SBATCH --gres=gpu:4
#SBATCH -p debug
#SBATCH --qos=restrained
#SBATCH --time=00:30:00
#SBATCH --exclusive
#SBATCH --no-requeue

module load cuda

# Source GROMACS (GMX) and PLUMED executables.
source /bgfs01/insite/utsab.shrestha/programs/gmx_plumed4/gromacs-2018.8/install_dir/bin/GMXRC.bash
source /bgfs01/insite/utsab.shrestha/programs/gmx_plumed4/plumed-2.5.4/sourceme.sh

#==============================================
#VARIABLES BELOW TO BE ADJUSTED/CHANGED BY USER
#==============================================
NumReplicas=20			#Total number of replicas to be used.

HremdMaxH=0.1                   #Maximum hours for production HREMD.
				#Keep "HremdMaxH" always slightly less than the walltime.
#==============================================


#=======================================
# VARIABLES BELOW USUALLY NEED NO CHANGE
#=======================================
RepEx=500			#Attempt replica exchange every 500 MD steps, i.e., 1 ps.

CheckPt=5			#Write checkpoint file for MD run every 5 minutes.
#=======================================


# Check to make sure that there are enough processes
if [[ $NumReplicas -gt $SLURM_NPROCS ]]; then
    echo "ERROR: Number of replicas greater than the number of tasks allocated for this job. This will disrupt MPI scheduling."
    exit 1
fi

touch plumed.dat
mpirun -np $NumReplicas -mca btl sm,self,tcp --oversubscribe gmx_mpi mdrun -cpi -cpt $CheckPt -maxh $HremdMaxH -plumed plumed.dat -replex $RepEx -hrex -multi $NumReplicas

rm *#*
