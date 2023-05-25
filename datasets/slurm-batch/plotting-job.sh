#!/bin/bash
#SBATCH -J plot			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=4         # cpus per node 
#SBATCH --partition=condo-cnattras         
#SBATCH --time=0-12:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=slurm-batch/slurm-err/plot.e%J	     
#SBATCH --output=slurm-batch/slurm-out/plot.o%J	     
#SBATCH --qos=condo
#SBATCH --mem=50G

module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
conda activate myhepenv

cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets
echo "Input directory: $1"
./PlotTFile $1

