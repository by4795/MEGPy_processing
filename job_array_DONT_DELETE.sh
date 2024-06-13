#!/bin/bash

#SBATCH -J      multiflux
#SBATCH -A ukaea-ap002-cpu
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --ntasks=56 ##56 ##4
##SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --mail-type=END,FAILED,BEGIN
##SBATCH --cpus-per-task=32
###SBATCH --array=1##--40
# #SBATCH --no-requeue


#SBATCH -p cclake
conda init bash
##mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment
source /home/ir-sidd1/.conda/envs/tglf #location of the virtual envuronment 
conda activate tglf

#application="python /home/ir-sidd1/rds/rds-ukaea-ap002-mOlK9qn0PlQ/ir-sidd1/TGLF/tglf_nn_gyrokit/pipeline_config.py"

#python /home/ir-sidd1/rds/rds-ukaea-ap002-mOlK9qn0PlQ/ir-sidd1/TGLF/tglf_nn_gyrokit/pipeline_config.py

python /home/ir-sidd1/rds/rds-ukaea-ap002-mOlK9qn0PlQ/ir-sidd1/TGLF/tglf_nn_gyrokit/aaron_config.py

JOBID=$SLURM_JOB_ID

cd "/home/ir-sidd1/rds/rds-ukaea-ap002-mOlK9qn0PlQ/ir-sidd1/TGLF/tglf_nn_gyrokit"

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

