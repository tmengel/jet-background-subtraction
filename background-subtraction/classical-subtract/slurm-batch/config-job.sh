#!/bin/bash
#SBATCH -J configMult			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=4          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-24:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/classical-subtract/slurm-batch/slurm-err/configMult.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/classical-subtract/slurm-batch/slurm-out/configMult.o%J	     
#SBATCH --qos=condo



module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/background-subtraction/classical-subtract
conda activate myhepenv
./Configure $1 $2 
