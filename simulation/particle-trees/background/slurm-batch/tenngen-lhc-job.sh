#!/bin/bash
#SBATCH -J TGLHC			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-06:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/background/slurm-batch/slurm-err/TGLHC.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/background/slurm-batch/slurm-out/TGLHC.o%J	     
#SBATCH --qos=condo

module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/background
conda activate myhepenv
hostname
echo $1
echo $2
echo $3
echo $4

./GenerateTennGenPbPb $1 $2 $3 $4
