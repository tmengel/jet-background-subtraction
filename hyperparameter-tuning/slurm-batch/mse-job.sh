#!/bin/bash
#SBATCH -J mse			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-24:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/hyperparameter-tuning/slurm-batch/slurm-err/mse.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/hyperparameter-tuning/slurm-batch/slurm-out/mse.o%J	     
#SBATCH --qos=condo


module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/hyperparameter-tuning
conda activate myhepenv

$if = $1 + 2
$ie = $2 + 2
$ig = $3 + 2 
export OUTFILE=/lustre/isaac/scratch/tmengel/jet-background-subtraction/hyperparameter-tuning/results/${SLURM_JOB_ID}.json
# echo job id
python3 mse_vs_parameters.py --hidden_nodes_start $1 $2 $3 --hidden_nodes_end $if $ie $ig --output $OUTFILE


