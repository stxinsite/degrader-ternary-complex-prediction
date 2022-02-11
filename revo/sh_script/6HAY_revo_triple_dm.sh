#!/bin/bash -x
#SBATCH --reservation=advsim
#SBATCH -p project                                                                                                                                                                                         
#SBATCH --qos=maxjobs                                                                                                                                                                                      
##SBATCH -o cIAP_tl7_amber_ff_runs/rebinding/slurm_logs/%x_%J.out                                                                                                                                          
##SBATCH -e cIAP_tl7_amber_ff_runs/rebinding/slurn_logs/%x_%J.err                                                                                                                                          
#SBATCH --nodes=1                                                                                                                                                                                          
#SBATCH --ntasks=1                                                                                                                                                                                         
#SBATCH --gres=gpu:4                                                                                                                                                                                       
#SBATCH --cpus-per-task=12                                                                                                                                                                                 
#SBATCH --mem=32gb                                                                                                                                                                                         
#SBATCH --time=168:00:00                                                                                                                                                                                   
#SBATCH --exclusive                                                                                                                                                                                        
#SBATCH -J production_rebinding_btk-protac__wepy__DYNAMITE                                                                                                                                                 #SBATCH --no-requeue 

source /bgfs01/common/we_envs/conda/bin/activate /bgfs01/common/we_envs/envs/wepy_newnodes

#module use /bgfs01/common/we_envs/modulefiles
#module load stx_wepy

HOME_DIR=$pwd


JOBNAME=${SLURM_JOB_ID}

SCRATCH_DIR=$HOME_DIR/'simulations'
WORK_DIR=$SCRATCH_DIR/$JOBNAME

mkdir $WORK_DIR

#LOG=$WORK_DIR/log

cp -r $HOME_DIR/6HAY_inputs/ $WORK_DIR
cp $HOME_DIR/python_scripts/6HAY_triple_distance_metric.py $WORK_DIR

cd $WORK_DIR

# ------------------------------                                                                                                                                                                            
# The code for this script                                                                                                                                                                                  
# ===============================================================================                                                                                                                           
echo "------------"# 1>> $LOG 2>> $LOG
echo "Running script"# 1>> $LOG 2>> $LOG
echo "===============================================================================" #1>> $LOG 2>> $LOG
                                                                                                                                                             
python -u 6HAY_triple_distance_metric.py 1 2000 10000 48 4

echo "===============================================================================" #1>> $LOG 2>> $LOG
echo "done with script"# 1>> $LOG 2>> $LOG
echo "------------"# 1>> $LOG 2>> $LOG
#echo ""   1>> $LOG 2>> $LOG
