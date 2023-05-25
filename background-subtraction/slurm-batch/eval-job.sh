#!/bin/bash
#SBATCH -J evaljob			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=2          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-48:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/slurm-batch/slurm-err/evaljob.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/slurm-batch/slurm-out/evaljob.o%J	     
#SBATCH --qos=condo
#SBATCH --mem=32G


module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction
conda activate myhepenv
./Evaluator $1 $2 $3