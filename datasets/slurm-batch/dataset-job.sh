#!/bin/bash
#SBATCH -J dataset			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=4         # cpus per node 
#SBATCH --partition=condo-cnattras         
#SBATCH --time=0-12:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=slurm-batch/slurm-err/dataset.e%J	     
#SBATCH --output=slurm-batch/slurm-out/dataset.o%J	     
#SBATCH --qos=condo
#SBATCH --mem=50G

module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
conda activate myhepenv

cd /lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets
echo "Input directory: $1"
echo "Output directory: $2"

./PrepareData -d $1 -o $2 -v 1

