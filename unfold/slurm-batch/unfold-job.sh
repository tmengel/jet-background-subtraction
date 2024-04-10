#!/bin/bash
#SBATCH -J unfoldJob			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-12:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold/slurm-batch/slurm-err/unfoldJob.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold/slurm-batch/slurm-out/unfoldJob.o%J	     
#SBATCH --qos=condo



module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/unfold
conda activate myhepenv
# hostname
ROOT_COMMAND="root -l -b -q 'Unfold.C(\"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\")'"
echo $ROOT_COMMAND
eval $ROOT_COMMAND