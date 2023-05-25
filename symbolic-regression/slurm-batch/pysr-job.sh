#!/bin/bash
#SBATCH -J pysr			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=4          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-24:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/symbolic-regression/slurm-batch/slurm-err/pysr.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/symbolic-regression/slurm-batch/slurm-out/pysr.o%J	     
#SBATCH --qos=condo
#SBATCH --mem=30G


module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/symbolic-regression
conda activate myhepenv
python3 fitPysr.py -f $1 -o $2 -m $3