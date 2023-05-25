#!/bin/bash
#SBATCH -J pythia			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-12:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/signal/slurm-batch/slurm-err/pythia.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/signal/slurm-batch/slurm-out/pythia.o%J	     
#SBATCH --qos=condo

module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/signal
conda activate myhepenv
hostname
echo $1
echo $2


./GeneratePythiaPP $1 $2 $3
