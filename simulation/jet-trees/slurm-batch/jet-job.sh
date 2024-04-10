#!/bin/bash
#SBATCH -J Jet			       #The name of the job
#SBATCH -A ACF-UTK0011               # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1         # cpus per node 
#SBATCH --partition=campus        
#SBATCH --time=0-01:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/jet-trees/slurm-batch/slurm-err/Jet.e%J	     
#SBATCH --output=/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/jet-trees/slurm-batch/slurm-out/Jet.o%J	     
#SBATCH --qos=campus

module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/jet-trees
conda activate myhepenv
echo "Submitting Jet Jobs."
./JetAnalysis $1 $2 $3 $4 $5 $6
echo "Job complete."
